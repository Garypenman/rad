# RAD: Reaction Analysis & Design

**A Declarative, Vectorized RDataFrame Framework for High-Performance Physics Analysis**

RAD is a header-only C++ framework designed for the Electron-Ion Collider (EIC) and other high-intensity experiments. It abstracts complex multi-particle combinatorics, kinematics, and Monte Carlo truth matching into a high-level declarative API, leveraging ROOT's **RDataFrame** for parallel execution.

### Key Features
* **Vectorized Core:** Replaces heavy objects (`TLorentzVector`) with Structure-of-Arrays (SoA) `ROOT::RVec<double>`, enabling compiler auto-vectorization and high cache locality.
* **Thread-Safe by Design:** Custom actions (`SnapshotCombi`, `Histogrammer`) are built for `ROOT::EnableImplicitMT()`, allowing lock-free parallel processing.
* **Declarative Syntax:** Define *what* particles you want, not *how* to loop over them.
* **JIT Compilation:** Leveraging Cling to compile selection logic at runtime for maximum performance.
* **Automated Matching:** Built-in handling for EDM4hep/Podio associations (Rec $\leftrightarrow$ Truth).

---

## ⚡ Quick Start

### Installation
RAD is header-only. Simply clone the repository and add it to your ROOT include path.

```bash
git clone [https://github.com/dglazier/rad.git](https://github.com/dglazier/rad.git)
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$(pwd)/rad/include
```

### Basic Usage Example
Here is a complete analysis steering macro (`Analysis.C`) for reconstructing $Y(4260) \to J/\psi(\to e^+e^-) \pi^+ \pi^-$ using ePIC simulation data.

```cpp
#include "ePICReaction.h"
#include "KinematicsProcessor.h"
#include "PhysicsSelection.h"
#include "Histogrammer.h"

void Analysis() {
    // 1. Enable Parallel Processing
    ROOT::EnableImplicitMT();

    // 2. Initialize Reaction & Data Source
    // ePICReaction handles Podio file reading and Rec-to-Truth matching automatically
    epic::ePICReaction rad_df{"events", "input_data/*.root"};
    rad_df.SetBeamsFromMC(0, 1); // Alias MC Beams to Rec namespace
    rad_df.SetupMatching();      // Traverse Podio Association Arrays

    // 3. Define Candidates (Name, MC_Role_ID, Filter, Columns)
    rad_df.SetParticleCandidates("ele", 6, rad::index::FilterIndices(11),  {"rec_true_pid"});
    rad_df.SetParticleCandidates("pos", 7, rad::index::FilterIndices(-11), {"rec_true_pid"});
    rad_df.SetParticleCandidates("pip", 4, rad::index::FilterIndices(211), {"rec_true_pid"});
    rad_df.SetParticleCandidates("pim", 5, rad::index::FilterIndices(-211),{"rec_true_pid"});

    // 4. Generate Combinations & Flag Signal
    rad_df.MakeCombinations();
    rad_df.DefineTrueMatchedCombi("rec_is_signal", rad::Rec());

    // 5. Configure Kinematics (The Processor)
    rad::KinematicsProcessor kine(&rad_df, rad::Rec());
    
    // Topology: Construct J/psi from e+e-, then Y from J/psi pi+pi-
    kine.Creator().Sum("Jpsi", {{"ele", "pos"}});
    kine.Creator().Sum("Y",    {{"Jpsi", "pip", "pim"}});
    
    // Physics Variables
    kine.Mass("MassJ", {"Jpsi"});
    kine.Mass("MassY", {"Y"});
    
    // Initialize the Chain
    kine.Init();

    // 6. Physics Selection (Lazy Masking)
    rad::PhysicsSelection sel(kine);
    sel.AddCutRange("Y_Signal_Region", "MassY", 3.8, 4.6);
    sel.Compile();

    // 7. Save Flat Tree (SnapshotCombi)
    // "rec_is_signal" allows you to filter pure signal in the output
    std::vector<std::string> vars = {
        "rec_MassY", "rec_MassJ", 
        "rec_ele_px", "rec_ele_py", "rec_ele_pz", // Auto-generated component vectors
        "rec_is_signal" 
    };
    
    rad_df.BookSnapshotCombi("Output.root", "tree", vars, sel.GetMaskColumn());
    rad_df.TriggerSnapshots();
}
```

---

## 🏗 Architecture

### 1. The Kinematics Processor
Unlike traditional tools that use lists of objects (`std::vector<Particle>`), RAD uses **Structure of Arrays**. The `KinematicsProcessor` is a functor that consumes a map of indices (`RVecIndices`) and operates on raw data columns (`Px`, `Py`, `Pz`).

* **Zero-Copy Logic:** Combinatorics are just integer lookups.
* **Vectorization:** Operations like `Mass = sqrt(E^2 - P^2)` are applied to the entire column at once, utilizing CPU SIMD instructions.
* **Master/Clone:** You can easily process **Truth** (Monte Carlo) data in parallel by cloning the processor:
    ```cpp
    auto kine_truth = kine_rec.CloneForType("tru"); // Automatically maps "rec_ele" -> "tru_ele"
    kine_truth->Init();
    ```

### 2. Parallel I/O: `SnapshotCombi`
Writing TTrees from multiple threads usually causes segmentation faults. RAD solves this with **`SnapshotCombi`**:
* **Lock-Free:** Pre-allocates thread-local buffers and `TTree` instances based on the thread pool size.
* **Context Safety:** Guards `gDirectory` to prevent file operations from crashing the input reader.
* **Lazy Flattening:** Automatically flattens the nested combinatorial structure (`Event` $\to$ `Combination`) into a flat N-Tuple (one row per candidate).

### 3. Physics Selection (Masking)
Standard `Filter()` drops events, which is bad for sideband analysis. RAD uses **Masking**:
* `PhysicsSelection` creates a boolean column (mask).
* Histograms and Trees only fill entries where `mask == true`.
* This allows you to run multiple selections (Signal vs Sideband) in a single pass.

---

## 🧩 ePIC & Podio Integration

RAD provides specialized support for the EIC's **ePIC** detector and the **EDM4hep** data model.

### Truth Matching
In Podio, Rec-to-Truth links are stored in separate association arrays (`ReconstructedParticle2MCParticle`).
* **`ePICReaction`** automatically traverses these arrays using `rad::PodioAssociations`.
* It maps the unique Object IDs to finding the correct MC particle index.
* Use `SetupMatching()` to enable this. It allows you to define "Roles" (e.g., "The scattered electron is MC Particle #2") and flag reconstructed candidates that match these roles (`DefineTrueMatchedCombi`).

### Auxiliary Data (Detector Info)
You can propagate low-level detector info (Cluster IDs, Track Chi2) through the combinatorial engine using **Auxiliary Variables**.

```cpp
// 1. Define the raw column (e.g., Calorimeter Energy)
rad_df.DefineAssociation("clusters", {"EcalBarrelClusters"}, "energy");
rad_df.DefineProjection("rec_cal_energy", "rec_clusters_energy", "rad::util::First");

// 2. Pass it to the Processor
// This ensures "rec_ele_cal_energy" is created and aligned with the electron candidates
kine.PreModifier().SetMomentumFrom("ele", "rec_cal_energy");
```

---

## 📁 Output Structure

The output `TTree` is **flat** and analysis-friendly.
* **Scalars:** Event-level variables (e.g., Truth Beam Energy) are broadcasted (repeated) for every combination.
* **Vectors:** Combinatorial variables are flattened.

| Entry | rec_MassY | rec_ele_px | rec_is_signal | tru_BeamE |
| :--- | :--- | :--- | :--- | :--- |
| 0 | 4.21 | 1.5 | 1 | 18.0 |
| 1 | 3.85 | 1.2 | 0 | 18.0 |
| 2 | 4.10 | 1.8 | 0 | 18.0 |

*(Note: "Entry" here represents one candidate combination. In the raw RDataFrame, these are elements of an RVec)*.
# RAD Project Style Guide

This document defines the coding standards and architectural patterns for the **RAD (Reaction Analysis & Design)** framework. Adherence to this guide ensures thread safety, high performance with ROOT RDataFrame, and maintainability.

---

## 1. File Structure & Organization

### 1.1 Header-Only Implementation
* All code is contained within `.h` files.
* **Strict Separation of Interface and Implementation:**
    * **Class Body:** Must contain only member variable and method declarations.
    * **Implementations:** Must be placed outside the class body, at the bottom of the file.
* **Inline Keyword:** All implementations must be marked `inline` to prevent One Definition Rule (ODR) violations when included in multiple translation units.

**Example:**

    namespace rad {
        class MyClass {
        public:
            // Declaration only
            void DoSomething(int value);
        private:
            int _internalValue;
        };

        // Implementation at bottom
        inline void MyClass::DoSomething(int value) {
            _internalValue = value;
        }
    }

### 1.2 Headers and Namespaces
* All code must reside within the `rad` namespace (or sub-namespaces like `rad::io`, `rad::physics`).
* Use `#pragma once` for header guards.
* **Minimize Includes:** Do not include headers inside class headers unless types are explicitly required in declarations. Forward declare if possible.
* **Never** utilize `using namespace` in a header file (global scope pollution).

### 1.3 Common Types & Definitions
* **Prefer Standard Aliases:** Use the type aliases defined in `CommonDefines.h` rather than raw ROOT/STL types. This improves readability and facilitates global precision changes (e.g., swapping `Double_t` for `Float_t`).

* **Indices:**
    * Use `Indices_t` instead of `ROOT::RVecI` or `std::vector<int>`.
    * Use `RVecIndices` (or `RVecIndexMap`) for nested structures (`RVec<RVecI>`).

* **Kinematics Data:**
    * Use `RVecResultType` instead of `ROOT::RVec<double>` or `RVecD` for kinematic columns.
    * Use `ResultType_t` instead of `double` for scalar calculations.

* **Strings:**
    * Use `ParticleNames_t` instead of `std::vector<std::string>`.

**Example:**

    // Good
    void Calculate(const RVecIndices& map, const RVecResultType& px);

    // Avoid
    void Calculate(const ROOT::RVec<ROOT::RVecI>& map, const ROOT::RVec<double>& px);

---

## 2. Naming Conventions

* **Namespaces:** `lower_case` (e.g., `rad`, `rad::io`).
* **Classes:** `PascalCase` (e.g., `KinematicsProcessor`, `SnapshotCombi`).
* **Methods:** `PascalCase` (e.g., `InitMap`, `GetReactionIndex`, `FinalizeTask`).
* **Member Variables:** `_camelCase` with leading underscore (e.g., `_fileName`, `_totalCount`, `_merger`).
* **Local Variables:** `snake_case` or `camelCase` (be consistent within scope).
* **Template Parameters:** `T` or `PascalCase` (e.g., `Tp`, `RecordType`).
* **Types/Aliases:** `PascalCase` usually ending with `_t` (e.g., `IndexMap_t`, `ParticleNames_t`).
* **Constants/Enums:** `PascalCase` or `kPascalCase` (e.g., `ColType::Double`, `OrderX()`).

---

## 3. Documentation (Doxygen)

### 3.1 Requirement
All public methods and classes must have full Doxygen documentation using the Javadoc style `/** ... */` block before the declaration.

### 3.2 Required Tags
* **`@brief`**: A single-line summary.
* **`@details`**: (Recommended) Expanded explanation of logic, especially for complex RDataFrame operations or threading models.
* **`@param`**: Description of every argument.
* **`@return`**: Description of the return value (if not void).
* **`@throw`**: (Optional) Document specific exceptions (e.g., `std::runtime_error`) if the method performs critical checks.

**Example:**

    /** * @brief Forces an input particle to be registered in the map.
     * @details Essential for Master processors when a Linked Clone needs a particle 
     * that the Master itself does not explicitly use.
     * @param name The name of the input particle to require.
     * @return True if registration was successful.
     */
    bool RequireParticle(const std::string& name);

---

## 4. RDataFrame & Architecture Patterns

### 4.1 Container Standardization
* **Prefer `ROOT::VecOps::RVec<T>`** (or the aliases in `CommonDefines.h`) over `std::vector<T>` for all data that interacts with RDataFrame columns.
* This ensures seamless interoperability with JIT-compiled actions and vector arithmetic.

    // Good
    void Process(const RVecResultType& px);

    // Avoid (unless internal logic only)
    void Process(const std::vector<double>& px);

### 4.2 Thread Safety & ImplicitMT
* **Assume Multi-threading:** All Actions and Processors must be thread-safe by design to support `ROOT::EnableImplicitMT()`.
* **Lock-Free Design:** Avoid `std::mutex` inside `Exec()` or `operator()`. Use **Thread-Local Storage** or **Pre-allocated Vectors** indexed by `slot` (obtained via `InitTask`).
* **Pre-allocation Pattern:**

    // In Initialize():
    const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
    _threadLocalData.resize(nSlots); // Pre-allocate to avoid race conditions

### 4.3 File I/O & Resources
* **TDirectory Context:** When opening ROOT files in worker threads (e.g., `InitTask`), **always** guard the global state to prevent crashing the Input Reader.

    void InitTask(TTreeReader*, unsigned int slot) {
        TDirectory::TContext c; // Scoped guard to prevent gDirectory pollution
        _files[slot] = _merger->GetFile();
    }

* **Memory Stability:** Use `std::unique_ptr` for buffers linked to `TTree::Branch`. This prevents memory address invalidation if the parent container resizes.

    std::vector<std::unique_ptr<Buffer_t>> _buffers;

### 4.4 Lifecycle Management
Custom RDataFrame Actions (`RActionImpl`) must strictly follow the ROOT lifecycle to ensure correct initialization and data flushing:
1.  **Constructor:** Configuration only. No resource allocation.
2.  **`Initialize()`:** Open files, allocate mergers, pre-size vectors (Main Thread).
3.  **`InitTask(slot)`:** Connect thread to resources / Lazy initialization (Worker Thread).
4.  **`Exec(slot, ...)`:** High-performance event loop (Worker Thread).
5.  **`FinalizeTask(slot)`:** Flush thread-local buffers to merger (Worker Thread).
6.  **`Finalize()`:** Merge results, close files, and cleanup (Main Thread).

---

## 5. Kinematics & Physics Logic

### 5.1 Vectorization (Structure of Arrays)
* Do not create arrays of objects (e.g., `vector<TLorentzVector>`).
* Use **Structure of Arrays (SoA)**: separate `RVecResultType` for Px, Py, Pz, Mass.
* This enables CPU auto-vectorization and better cache locality.

### 5.2 Optimization
* **Avoid Redundant Math:** Do not compute Energy ($E = \sqrt{P^2 + M^2}$) inside combinatorial loops. Store and copy **Mass** instead. Compute Energy only when plotting/saving if strictly necessary.
* **Functors:** Implement processors as stateless functors (`operator()`) that consume `RVecIndices` rather than objects.

---

## 6. Error Handling & Safety

### 6.1 Map Access
* Never use raw `map::at()` or `[]` operators on critical maps (like particle indices) without a wrapper check.
* Use a "Safe Getter" pattern that throws a `std::runtime_error` with a descriptive message including context (Prefix, Suffix, missing key).

### 6.2 Empty File Recovery
* Output actions must handle the case where **0 events** pass selections.
* Always ensure the output ROOT file contains the expected `TTree` header (even if empty) to prevent downstream crashes.

### 6.3 Explicit Triggering
* Because RDataFrame is lazy, complex chains (especially Snapshots) should be explicitly triggered via `rad_df.TriggerSnapshots()` before the macro exits to guarantee data is flushed.

---

## 7. Framework Design Patterns

### 7.1 Master/Clone Pattern
* **Master:** Configures the primary topology and definitions.
* **CloneForType:** Copies configuration to a different data stream (e.g., Rec -> Truth) by swapping prefixes.
* **CloneLinked:** Creates a new hypothesis on the same data stream by adding a unique suffix.

### 7.2 Group Porting
* When cloning across types (Rec -> Truth), use `PortGroupsFrom` to inspect the parent's groups and attempt to register equivalent groups for the new type automatically.

### 7.3 Redefinition
* Use `Redefine()` when passing a C++ lambda or function pointer.
* Use `RedefineExpr()` when passing a string containing C++ code (JIT compilation).

### 7.4 Lambda Type Safety
* Prefer explicit argument types in lambdas over `auto` when dealing with `RVec` to ensure ROOT can infer types correctly during compilation.
* **Example:** `[](const Indices_t& a) { ... }` instead of `[](auto a) { ... }`.
---
applyTo: '**'
---

# DEVELOPING GUIDELINES

1. **NO template, NO lambda, No PCH, No macro define (Use constexpr if possible for type safety)**
2. **Plain design**, no extra namespaces more than `hpcfem` (To seperate our own lib from the third-party ones).
3. **Always check to avoid enum conflicts.**
4. **No code nesting.** No struct/class/enum inside struct/class.
5. **Naming uniqueness.** Everytime you introduce new class/enum/struct, document them in `docs/hpcfem-doc`, and check to avoid duplication before adding new ones.
6. **Use Doxygen-style comments.** The comments must be clean, with only interface description, functionalities and TODOs(If any). **NEVER MARK PROGRESS IN COMMENTS.**
7. Guarantee that all the examples/tests/benchmarks can all compile and run smoothly before a `commit`
8. If for anything you only made interface or simplified implementation, please EXPLICITLY mark it as `TODO` in comment. And you should clear up all the TODOs before a `merge`.
9. **use `git` for version management, usually you work on a new branch with name related to your work, after a period of work is done, merge it into `dev` and delete it.** If a large scale refactor is needed, please open an `unstable` branch to do it. After refactor completes, merge it to the `dev` branch. NEVER TOUCH THE `master` or `main` BRANCH!
10. **Adhere to TDD principle.** Recommended workflow: PLANNING->WRITING TESTS->IMPLEMENTING(for test passing)->TESTING->REFACTORING(for better code quality)->DOCUMENTING->STRESS TEST(For capability with more massive data scale, e.g. 1e7 Dofs / 1e8 element mesh). 
11. **What if you get stuck:** If you'd try more than 3 times but still fail on the same problem, you MUST stop and try to figure out the root cause from a higher perspective. You must have a project-level view instead of struggling repetitively for a in-place fix.
12. **What kind of implement to choose:** When you have multiple routines to implement the same functionalities, choose by the following principles: Testablility>Scalability>Performance>Maintainability>Consistency>Simplicity.
13. **YAGNI(You Ain't Gonna Need It):** While generality is important, you should never introduce unnecessary abstraction and complexity, unless the goal itself requires so.
14. **Keep project code readable and maintainable:** No more than 500 lines in a single file. If you have to do so, it usually prove to be a flaw in design and you may have to do a refactor.
15. **NO MAGIC NUMBER OR CONSTANT:** We allow no magic number or other kind of constants that comes from nowhere. All the constants must be defined with meaningful names at the begining of the file or in a separate config file.
16. **Naming Traditions**:
   - Class/Struct/Enum: PascalCase
   - Function/Method: camelCase
   - Variable/Member: camelCase
   - Constant: UPPER_SNAKE_CASE
   - File/Folder: lowercase_with_underscores/Keep the same name as the main class inside
   - Test Case: test_case_description_with_underscores
17. **Put .hpp files and .cpp files in the same location in the `src`**
18. **Clean include path:** No use of relative path in include like `../`, instead, always configure the include path property in CMakeLists.txt
19. **Correct Python configuration if needed:** Use the conda environment `hpc-fem-playground` for python use.
20. **Proper problems setup:** All the problem setup in tests and examples should be aligned to pratical use (3D cases preferred, vector field preferred).
21. **Build correctly:** Debug build in `cmake-build-debug`, Release build in `cmake-build-release`.
22. **Plan and checklists:** For the future developers to take on your work, please always provide a plan document with clear todos and checklists. Before implementing anything, check out the document file. Whenever you finished something in the checklist, please update the document.
23. **Proper cmake configuration:** Use CPM to mangage dependencies. After any update on cmake scripts, you should rerun configuration before rebuild.
24. **Properly arrage the file:** Traditionally you shouldn't put more than 15 files in a directory. If you're doing so, you should consider refactoring.
25. **Consider MPI parallelism:** Your code should be able to run in parallel with MPI without breaking down. Everytime you add a new functionality, your test shall include both the serial and parallel cases.
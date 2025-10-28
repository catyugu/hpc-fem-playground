---
applyTo: '**'
---

# DEVELOPING GUIDELINES

1. **NO template, NO lambda, No PCH, No macro define (Use constexpr if possible for type safety)**
2. **Plain design**, no extra namespaces more than `hpcfem` (To seperate our own lib from the third-party ones).
3. **Always check to avoid enum conflicts.**
4. **No code nesting.** No struct/class/enum inside struct/class.
5. **Naming uniqueness.** Everytime you introduce new class/enum/struct, document them in `docs/hpcfem-doc`, and check to avoid duplication before adding new ones.
6. **Use Doxygen-style comments.**
7. Guarantee that all the examples/tests/benchmarks can all compile and run smoothly before a `commit`
8. If for anything you only made interface or simplified implementation, please EXPLICITLY mark it as `TODO` in comment. And you should clear up all the TODOs before a `merge`.
9. **use `git` for version management, usually you work on `dev` branch.** If a large scale refactor is needed, please open an `unstable` branch to do it. After refactor completes, merge it to the `dev` branch. For every markstone you've reached, please merge from `dev` to `master`.
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
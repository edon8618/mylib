# LinkCutTree verification

Run from repository root. No network or other dependencies are needed for the main tests.

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic tests/link_cut_tree/test.cpp -o /tmp/lct_test
/tmp/lct_test 100
g++ -std=c++20 -O2 tests/link_cut_tree/api_compile.cpp -o /tmp/lct_api
/tmp/lct_api
g++ -std=c++20 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer tests/link_cut_tree/test.cpp -o /tmp/lct_san
ASAN_OPTIONS=detect_leaks=0 /tmp/lct_san 10
python3 tests/link_cut_tree/run_bench.py tests/link_cut_tree/benchmark_200000.csv 200000 3
python3 tests/link_cut_tree/run_bench.py tests/link_cut_tree/benchmark_500000.csv 500000 3
python3 tests/link_cut_tree/summarize_bench.py
```

`test.cpp` compares all requested APIs with BFS on a naive forest, using signed integer
addition and affine composition modulo 998244353 (including zero multipliers).
`api_compile.cpp` verifies the absence of 3-argument endpoint overloads for a
distinguishable value type, and supports a monoid with no default constructor.
For S=int, the 3-integer ID overload cannot be distinguished from accidental endpoints.
The seed count controls 10000 iterations per seed per monoid. Assertions must remain enabled.
It also checks opening-direction semantics, deletion/reinsertion, identity-edge materialization,
legacy APIs, and a 200000-vertex chain. Node reuse cannot recycle edge IDs.

ASan and UBSan are enabled together. Leak detection is disabled because LeakSanitizer cannot
run under the execution environment's ptrace configuration; this is not a leak-check result.

## Benchmarks

`old_LinkCutTree.hpp` is the exact original header from commit
`c64718ccd2586e5391177c174c448c0067be566e` of edon8618/mylib.
It is retained solely as a reproducible reference and depends on the benchmark's standard
headers and `using namespace std`, as in the original repository.

Compiler: g++ -std=c++20 -O2 -DNDEBUG, with -DOLD and/or -DAFFINE.
No architecture-specific flags. The program's fixed PRNG seed is 1234567.
N=Q, with 25% updates and 75% random endpoint path products. Modes:

| Mode | Workload |
|---|---|
| 0 | Vertex values, random recursive tree |
| 1 | Symmetric edge values, random recursive tree |
| 2 | Directional edge values, random recursive tree |
| 3 | Vertex and directional edge values, random recursive tree |
| 4 | Vertex values, long path |
| 5 | Vertex values, random recursive tree; 25% random-edge cut/relink pairs |

Mode 5 reconnects the cut edge's original endpoints, so the topology stays a tree.
It exercises dynamic operations but does not represent all possible topology changes.
Random correctness tests do change topology and split/merge components.

The old vertex workloads use their native vertex-only representation. For symmetric edge
values, the old version manually allocates N-1 additional vertices. Directional edge modes
have no old baseline because the original API cannot encode directional values directly.
Old uses split/get_prod, new uses prod; RNG streams and comparable checksums must agree.
Both affine and integer-sum monoids run in all modes. Integers are `long long` in benchmarks,
`int` in correctness tests. Each binary is a separate process; new/old ordering alternates.
CSV records build time and operation time separately (wall time and process CPU time),
peak process RSS from getrusage, and checksum. RSS includes the program and its input arrays.
Use medians over 3 runs. CPU times help distinguish scheduling delays from computation;
results are specific to this shared execution environment, not a judge time guarantee.

## Currency network example

`examples/currency_network.cpp` uses ACL modint and documents an assumed input format,
because no official input specification was supplied. It is a runnable API example,
not a verified submission to the unspecified original problem.

```sh
g++ -std=c++20 -O2 -I/path/to/ac-library examples/currency_network.cpp -o /tmp/currency
/tmp/currency < tests/link_cut_tree/currency.in > /tmp/currency.out
diff -u tests/link_cut_tree/currency.out /tmp/currency.out
```

The fixed fixture verifies disconnected queries, both directions, reversed update endpoints,
opening-direction updates, deletion, and reversed reinsertion. The six documentation examples
are also compiled and executed as part of the reported validation.

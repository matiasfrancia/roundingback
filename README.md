# RoundingBack

**RoundingBack** is a native backbone extractor for Pseudo-Boolean Optimization (PBO) instances, built upon the internals of [RoundingSat]([https://github.com/roundingsat/roundingsat](https://gitlab.com/MIAOresearch/software/roundingsat/)).
It reuses the solver's native data structures and algorithms to perform backbone extraction at low level.

This tool was developed for the experiments presented in a research paper submitted to the **SAT 2025 Conference**.

---

## Compilation

To build the backbone extractor, run the following commands in the root directory:

```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

For a debug build:

```bash
cd build_debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

RoundingBack supports optional integration with the SoPlex LP solver via RoundingSat.
To enable it, download soplex and build with:

```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Dsoplex=ON ..
make
```

For full instructions on how to build with Soplex, refer to the [RoundingSat](https://gitlab.com/MIAOresearch/software/roundingsat/) repository.

---

## Usage

Once the Release version of the backbone extractor is built, to extract the backbone of a PBO instance run:

```bash
./extract_backbone.sh <path_to_instance.opb>
```

If the instance is an optimization problem, you can optionally pass the optimum value to skip solving it:

```bash
./extract_backbone.sh <path_to_instance.opb> --opt <optimum_value>
```

You can also set a total timeout (in seconds), default is 3600 seconds:

```bash
./extract_backbone.sh <path_to_instance.opb> --timeout <seconds>
```

Both `--opt` and `--timeout` can be combined:

```bash
./extract_backbone.sh <path_to_instance.opb> --opt <optimum_value> --timeout 3600
```

If no optimum is passed, the script will attempt to detect whether the instance is decisional or optimization, extract the optimum value if needed, and compute the backbone accordingly, all within the total time budget.

---

## Citation

RoundingBack builds upon the core components of **RoundingSat**, a state-of-the-art solver for 0-1 Integer Linear Programs. If you use RoundingBack or RoundingSat in your work, please consider citing the following papers:

- **[EN18]** J. Elffers, J. Nordström. *Divide and Conquer: Towards Faster Pseudo-Boolean Solving*. In Proceedings of IJCAI 2018. 

- **[DGN20]** J. Devriendt, A. Gleixner, J. Nordström. *Learn to Relax: Integrating 0-1 Integer Linear Programming with Pseudo-Boolean Conflict-Driven Search*. In Proceedings of CPAIOR 2020 / Constraints Journal.  

- **[D20]** J. Devriendt. *Watched Propagation for 0-1 Integer Linear Constraints*. In Proceedings of CP 2020.  

- **[DGDNS21]** J. Devriendt, S. Gocht, E. Demirović, J. Nordström, P. J. Stuckey. *Cutting to the Core of Pseudo-Boolean Optimization: Combining Core-Guided Search with Cutting Planes Reasoning*. In Proceedings of AAAI 2021.  

---

## License

This project reuses and extends the codebase of RoundingSat, which is licensed under the **MIT License**. See the LICENSE file for details. Contributions by Matías Francia are also released under the MIT License.

---

## Related Dataset

We provide a dataset of backbone extraction results from this tool at:
[https://github.com/matiasfrancia/dataset-pbo-backbones](https://github.com/matiasfrancia/dataset-pbo-backbones)

It includes over 8000 extractions on instances from the OPT-LIN track of the Pseudo-Boolean Competition 2024.


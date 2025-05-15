# Quantum Hartley Transform (QHT) Implementation

This repository provides a Qiskit implementation of the Quantum Hartley Transform (QHT) algorithm as described in the paper ["Public-Key Quantum Money and Fast Real Transforms"](https://arxiv.org/abs/2503.18890). The implementation follows the pseudocode in Section 4.

---

## Correspondence to Paper

| Code Module / File   | Paper Section            | Algorithm           | Page  | Notes                                    |
| -------------------- | ------------------------ | ------------------- | ----- | ---------------------------------------- |
| `qht_recursion.py`   | Section 4                | Algorithm 1 (𝒬𝐻𝒯)  | 7–10 | Translation of the recursive pseudocode  |

---

## Command-Line Interface

    # Default behavior (n=3)
    python qht_recursive.py

## Installation

Install Python dependencies:

    pip install -r requirements.txt

## Circuit Diagram for n=3

![Quantum Hartley Transform Circuit](img/circuit.png)

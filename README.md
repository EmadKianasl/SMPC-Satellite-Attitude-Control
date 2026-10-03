# Privacy-Preserving Distributed Satellite Attitude Control

## Overview

This project investigates **privacy-preserving distributed control for satellite attitude systems** using **Secure Multiparty Computation (SMPC)** techniques.

The main idea is to perform attitude control while protecting sensitive sensor and control information from individual computing parties. The proposed framework combines **quantized sensor data**, **finite-field arithmetic**, **encrypted signal processing**, and **Shamir's Secret Sharing** to enable secure computation without directly exposing the underlying private data.

The implementation is developed using **Python/Jupyter, SageMath, and MATLAB**.

---

## Key Features

* 🔐 Privacy-preserving satellite attitude control
* 🤝 Secure Multiparty Computation (SMPC)
* 🔢 Quantization of sensor measurements
* 🧮 Finite-field arithmetic for secure computation
* 🔒 Privacy-preserving A/D conversion
* 🔓 Privacy-preserving D/A conversion
* 🧩 Shamir's Secret Sharing
* 🛰️ Satellite attitude-control application
* 🐍 Python / Jupyter implementation
* 📐 SageMath implementation
* 📊 MATLAB-based modeling and analysis

---

## Research Motivation

Modern spacecraft increasingly rely on distributed computing, cloud-based processing, and networked control architectures. While these approaches can provide computational and architectural benefits, they can also introduce privacy and information-security concerns.

In a conventional control architecture, sensor measurements and controller signals may be directly available to the computing platform. This can expose sensitive information about the spacecraft state or control actions.

This project explores an alternative architecture in which sensitive information is protected during computation.

The objective is to investigate how cryptographic techniques can be integrated with control-system computation while preserving the functionality of the satellite attitude controller.

---

## Conceptual Architecture

```text
              Satellite Sensors
                     │
                     ▼
             ┌────────────────┐
             │ Quantization   │
             └───────┬────────┘
                     │
                     ▼
             ┌────────────────┐
             │ Secure A/D     │
             │ Conversion     │
             └───────┬────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │ Secret Sharing /     │
          │ Secure Encoding      │
          └──────────┬───────────┘
                     │
          ┌──────────▼───────────┐
          │ Secure Multiparty    │
          │ Computation (SMPC)   │
          └──────────┬───────────┘
                     │
                     ▼
             ┌────────────────┐
             │ Attitude       │
             │ Controller     │
             └───────┬────────┘
                     │
                     ▼
          ┌──────────────────────┐
          │ Secure D/A           │
          │ Conversion           │
          └──────────┬───────────┘
                     │
                     ▼
             Satellite Actuation
```

The architecture separates the acquisition, secure representation, distributed computation, and control stages so that sensitive information does not need to be exposed directly to an individual computational party.

---

## Secure Multiparty Computation

The project applies **Secure Multiparty Computation (SMPC)** to the control problem.

In an SMPC framework, multiple parties can collaboratively perform computations on protected information without requiring every party to access the original private values.

For this project, the secure computation framework is based on:

* Secret sharing
* Finite-field arithmetic
* Quantized signals
* Secure distributed computation
* Reconstruction of computed results

This provides a bridge between **control engineering** and **cryptographic computation**.

---

## Shamir's Secret Sharing

A central component of the privacy-preserving architecture is **Shamir's Secret Sharing**.

Instead of directly transmitting a sensitive value, the value can be represented through shares distributed among multiple parties.

Conceptually:

```text
                 Private Value
                      │
                      ▼
             ┌─────────────────┐
             │ Shamir Secret   │
             │ Sharing         │
             └────────┬────────┘
                      │
        ┌─────────────┼─────────────┐
        ▼             ▼             ▼
     Share 1       Share 2       Share 3
        │             │             │
        ▼             ▼             ▼
     Party 1        Party 2        Party 3
        │             │             │
        └─────────────┼─────────────┘
                      ▼
              Secure Computation
                      │
                      ▼
               Result Reconstruction
```

The use of secret sharing allows the distributed parties to participate in the computation without requiring the original sensitive value to be explicitly revealed to each party.

---

## Quantized Sensor Data

The framework considers **quantized sensor measurements** as part of the privacy-preserving control pipeline.

Quantization converts continuous-valued measurements into a discrete representation suitable for subsequent processing and finite-field operations.

A conceptual quantizer can be represented as

$$
q(x) = Q(x)
$$

where:

* \(x\) is the original sensor measurement,
* \(Q(\cdot)\) is the quantization operation,
* \(q(x)\) is the resulting discrete representation.

The quantized representation can then be mapped into an appropriate finite-field representation for secure computation.

---

## Finite-Field Arithmetic

Finite-field arithmetic is used to support secure numerical operations.

Instead of performing the required operations directly over ordinary real-valued numbers, data can be represented within a finite field:

$$
\mathbb{F}_p
$$

where \(p\) is an appropriate prime modulus.

Operations are therefore performed modulo \(p\):

$$
(a+b)\bmod p
$$

and

$$
(a b)\bmod p.
$$

This representation is particularly useful for cryptographic protocols and secret-sharing-based computation.

---

## Privacy-Preserving A/D and D/A Conversion

An important aspect of the project is the treatment of the **Analog-to-Digital (A/D)** and **Digital-to-Analog (D/A)** interfaces within the privacy-preserving architecture.

The project investigates finite-field-based processing for:

### A/D Conversion

```text
Physical Sensor
      │
      ▼
Measurement
      │
      ▼
Quantization
      │
      ▼
Finite-Field Representation
      │
      ▼
Secure Processing
```

### D/A Conversion

```text
Secure Controller Output
          │
          ▼
Finite-Field Representation
          │
          ▼
Reconstruction
          │
          ▼
D/A Conversion
          │
          ▼
Actuator
```

This allows the secure computation framework to be connected to the conventional physical control system.

---

## Satellite Attitude Control

The application considered in this project is **satellite attitude control**.

Satellite attitude control is responsible for regulating the spacecraft's orientation and generating appropriate control actions based on attitude measurements.

The project focuses specifically on integrating privacy-preserving computation into this control architecture.

The overall concept can be summarized as:

$$
\text{Sensor Measurements}
\rightarrow
\text{Secure Representation}
\rightarrow
\text{Distributed Computation}
\rightarrow
\text{Control Command}
\rightarrow
\text{Satellite}
$$

This creates a connection between:

* Aerospace control
* Distributed control
* Cryptography
* Secure computation
* Digital control implementation

---

## Software and Tools

The project uses multiple computational environments:

| Tool                 | Role                                                     |
| -------------------- | -------------------------------------------------------- |
| **Python**           | Numerical and algorithmic implementation                 |
| **Jupyter Notebook** | Interactive development and experimentation              |
| **SageMath**         | Finite-field and mathematical/cryptographic computations |
| **MATLAB**           | Control-system modeling, simulation, and analysis        |

---

## Repository Structure

The current repository contains the main project implementation under:

```text
SMPC-Satellite-Attitude-Control/
│
├── SMPC Satellite’s Attitude Control/
│
└── README.md
```

The implementation folder contains the project materials associated with the Python/Jupyter, SageMath, and MATLAB portions of the work.

> The exact internal filenames are intentionally not listed here because the accessible GitHub repository view does not expose the contents of the nested project folder.

---

## Methodology

The overall methodology can be summarized in the following stages:

### 1. Sensor Measurement

Satellite attitude information is obtained from the sensing system.

### 2. Quantization

Continuous measurements are converted into discrete values suitable for digital and secure processing.

### 3. Finite-Field Representation

The quantized values are mapped into a finite-field representation.

### 4. Secret Sharing

Sensitive information can be distributed among multiple computational parties using **Shamir's Secret Sharing**.

### 5. Secure Computation

The parties perform the required computations using an SMPC framework without directly exposing the underlying private information.

### 6. Control Computation

The secure computational result is used to obtain the required attitude-control command.

### 7. Reconstruction and Actuation

The resulting control information is reconstructed and transferred toward the physical actuation stage.

---

## Research Topics

This project combines several research areas:

* **Satellite Attitude Control**
* **Distributed Control Systems**
* **Privacy-Preserving Control**
* **Secure Multiparty Computation**
* **Shamir's Secret Sharing**
* **Finite-Field Arithmetic**
* **Cryptographic Computation**
* **Quantized Control Signals**
* **Digital Control Systems**
* **Aerospace Engineering**

---

## Why This Project Is Interesting

The project explores an interdisciplinary problem at the intersection of **control engineering and cybersecurity**.

Traditional control-system design generally assumes that the controller can access the required measurements and perform computations directly. Privacy-preserving control introduces an additional requirement:

> **How can a control system perform useful computations while minimizing the exposure of sensitive information?**

The satellite-attitude-control application provides a concrete engineering scenario for studying this question.

---

## Getting Started

Clone the repository:

```bash
git clone https://github.com/EmadKianasl/SMPC-Satellite-Attitude-Control.git
cd SMPC-Satellite-Attitude-Control
```

Then open the project materials in the corresponding computational environment:

* **Jupyter / Python** for Python-based implementations
* **SageMath** for finite-field and secret-sharing computations
* **MATLAB** for control-system modeling and simulation

Because the repository contains implementations across multiple environments, the appropriate software should be selected according to the specific part of the project being studied.

---

## Project Scope

This repository is primarily an academic/research implementation demonstrating the integration of privacy-preserving computation with satellite attitude control.

The project should therefore be viewed as a research-oriented framework rather than a complete flight-qualified spacecraft control system.

Practical spacecraft deployment would require additional considerations such as:

* Real-time computational constraints
* Communication delays
* Numerical precision
* Fault tolerance
* Hardware implementation
* Cryptographic security assumptions
* Sensor and actuator limitations
* Space-qualified computing hardware

---

## Author

**Emad Kian Asl**

Mechanical Engineering / Control & Dynamics Projects

GitHub: [@EmadKianasl](https://github.com/EmadKianasl)

---

## License

Unless otherwise specified in the repository, please contact the author before using the project for commercial or redistributed applications.

---

## Acknowledgment

This project brings together concepts from **aerospace control, distributed computation, cryptography, and secure multiparty computation** to investigate privacy-preserving approaches to satellite attitude control.

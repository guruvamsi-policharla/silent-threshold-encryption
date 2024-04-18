# Silent Threshold Encryption [ePrint:2024/263](https://eprint.iacr.org/2024/263)

Rust implementation of the silent-threshold encryption introduced in [ePrint:2024/263](https://eprint.iacr.org/2024/263).

Use ```cargo bench``` to benchmark `setup`, `encryption`, and `decryption`.

Use ```cargo run --example endtoend``` to check correctness of the implementation.

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Overview
* [`src/setup`](src/setup.rs): Contains an implementation for sampling public key pairs and aggregating keys of a chosen committee. Also contains the `partial_decryption` method which is essentially a BLS signature. Note that the `get_pk` method runs in quadratic time. This can be reduced to linear time by preprocessing commitments to lagrange polynomials.
* [`src/encryption`](src/encryption.rs): Contains an implementation of the `encrypt` method for the silent threshold encryption scheme.
* [`src/decryption`](src/decryption.rs): Contains an implementation of `agg_dec` which gathers partial decryptions and recovers the message.

## License
This library is released under the MIT License.

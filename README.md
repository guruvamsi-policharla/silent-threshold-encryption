# Silent Threshold Encryption [ePrint:2024/263](https://eprint.iacr.org/2024/263)

Rust implementation of the silent-threshold encryption introduced in [ePrint:2024/263](https://eprint.iacr.org/2024/263). Benchmarks reported in the paper were run on a 2019 MacBook Pro with a 2.4 GHz Intel Core i9 processor. The library has been confirmed to work with version 1.76.0 of the Rust compiler. 

An end to end example is provided in the `examples/` directory.

## Dependencies
Install rust via:

```curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh```

## Benchmarking
Use ```cargo bench``` to benchmark `setup`, `encryption`, and `decryption`. This is expected to take approximately 20 minutes. To run a specific benchmark, use ```cargo bench --bench <bench_name>```.

Use ```cargo run --example endtoend``` to check correctness of the implementation.

The results are saved in the `target/criterion` directory. A concise HTML report is generated in `target/criterion/index.html` and can be viewed on a browser (Google Chrome recommended).

If you wish to benchmark for a different set of parameters, you can modify the files in the `benches/` directory. 

## Unit Tests
Additionally, you can find individual unit tests at the end of the respective files in the `src/` directory. These can be run using ```cargo test <test_name>```. This will allow you to test the correctness of the implementation.

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Overview
* [`src/setup`](src/setup.rs): Contains an implementation for sampling public key pairs and aggregating keys of a chosen committee. Also contains the `partial_decryption` method which is essentially a BLS signature. Note that the `get_pk` method runs in quadratic time. This can be reduced to linear time by preprocessing commitments to lagrange polynomials.
* [`src/encryption`](src/encryption.rs): Contains an implementation of the `encrypt` method for the silent threshold encryption scheme.
* [`src/decryption`](src/decryption.rs): Contains an implementation of `agg_dec` which gathers partial decryptions and recovers the message.

## License
This library is released under the MIT License.

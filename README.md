# secp256k1

Rust implementation of secp256k1 curve. Bare metal implementation without using any crates.  

# Running the project

```
cargo test
```

# Approach to tackle bigints

Using `Vec<u8>` vectors of hex for bigint replacement.

# Inspired from noble-secp256k1 implementation
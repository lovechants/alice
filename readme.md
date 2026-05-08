# Alice

A rigorous topology and algebra library from first principles, with Lie group machinery.
> on crates.io as `alice-math`. Import as `alice` via the package alias

```toml
[dependencies]
alice = { package = "alice-math", version = "0.1.0" }
```

## What it is

alice builds the mathematical foundations needed for Lie group simulation and Riemannian optimization from the ground up in Rust. The algebraic hierarchy starts at `Magma` and builds to `Field`. The topology layer encodes open set axioms, manifolds, and smooth structure. The groups module provides concrete matrix Lie groups with correct constraint enforcement and `no_std` compatibility.

## Modules

- **`alice::core`** — Algebraic hierarchy (`Magma` through `Field`), `FiniteF64`, `Rational`, morphisms
- **`alice::topology`** — Topological spaces, manifolds, smooth structure, metric spaces
- **`alice::algebra`** — Lie algebras, Lie bracket, Baker-Campbell-Hausdorff formula
- **`alice::maps`** — Exponential and logarithmic maps via Padé approximation and eigendecomposition
- **`alice::groups`** — `GL`, `SL`, `SO`, `SE`, `Aff` with full constraint enforcement at construction

## Quick example

```rust
use alice::groups::so::So;
use alice::maps::exp_log::HasExpMap;
use alice::topology::manifold::Dynamic;

let theta = std::f64::consts::PI / 4.0;
let rotation = So::new(
    vec![theta.cos(), -theta.sin(), theta.sin(), theta.cos()],
    Dynamic::new(2),
).unwrap();

let log = rotation.log().unwrap();
let recovered = So::exp(&log);
```

## Design

- `no_std + alloc` compatible
- Algebraic axioms enforced at the type level where possible, at construction otherwise
- All group constructors return `Result<Self, GroupError>` rather than panicking

## Companion crate

`riemannian-gd` provides Riemannian gradient descent over Lie groups, consuming `alice-math` as its foundation.

## Status

`0.1.0` — API is functional and tested. Atlas implementations for concrete groups, additional Lie group examples, and Manim visualizations are planned for subsequent releases.

## License

MIT

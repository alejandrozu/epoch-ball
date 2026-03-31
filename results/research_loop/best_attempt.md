# Strongest Current Path

- Score: 0.27
- Ready for paper: no
- Hypothesis: The repository's optimization-based polynomial Hamiltonian is not a correct general solution; it only supports a near-trivial epsilon regime and may not even outperform its own linear translate-and-scale seed.

## Why it is currently best

Installed CPU PyTorch so the provided code could be executed. Reproduced the shipped verifier pass for k=10, eps=0.99. Added src/hamiltonian_audit.py and used it to audit the shipped Hamiltonian on seeds 0, 123, and 999 with 1000 points per ball and 500 RK4 steps; the min sampled ratio was 0.012385 and the Hamiltonian call time stayed around 3-5 ms for a batch of 5000 points. Captured a timed trac

## Known limits

The optimization-based Hamiltonian family did not yield evidence for denser packing. The shipped optimized Hamiltonian underperformed the plain linear seed (min sampled ratio 0.012385 versus 0.013649 for the seed). The harder eps=0.98 run started at ratio 0.01357 and showed no best-ratio improvement through the recorded iterations. The included verifier is also non-independent because the generato

## Next move

Switch to a qualitatively different construction family: abandon global low-degree polynomial optimization and try a constructive localized packing method, such as toric/moment-polytope style embeddings or sequential compactly supported Hamiltonians, reusing src/hamiltonian_audit.py for verification.

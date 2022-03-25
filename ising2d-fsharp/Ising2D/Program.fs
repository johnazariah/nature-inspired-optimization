open Ising2D_Metropolis

let T = 0.9<K>
let num_iterations = 500_000
let (before, after, duration, flips) = Ising2D_Metropolis.Solve T num_iterations

printfn $"{before}"
printfn $"{after}"
printfn $"Solving Ising2D Metropolis with F#:\{L} x {L} ({num_iterations} iterations) at {T}K took {duration.TotalMilliseconds} ms. {flips} spins were flipped."
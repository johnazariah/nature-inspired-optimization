module Ising2D_Metropolis

[<Measure>] type K // temperature
[<Measure>] type E // hamiltonian energy

// It would be nice to make this a parameter to the module but the .NET type system doesn't have ways to define dependent types
let L = 20

// This is part of the 2D Ising specification to compute which cells interact as neighbours
//  Wrap around interactions happen with the modulo operator
let inline private neighbour (x, y, p) =
    match p with
    | 0 -> ((x + 1) % L), (y)
    | 1 -> ((x + L - 1) % L), (y)
    | 2 -> (x), ((y + 1) % L)
    | 3 -> (x), ((y + L - 1) % L)
    | _ -> failwith "how?"
    
// This is the problem definition for dimension LxL 
//    The JijMatrix represents a ferromagnet when all interactions are +/- 1
//    The JijMatrix represents an anti-ferromagnet when the interactions are alternately +1 and -1 to form a checkerboard
//    The JijMatrix represents a general spin-glass if the interations are otherwise varied and non homogenous
type private JijMatrix = JijMatrix of float<E>[,,]
with
    member inline this.Unapply = match this with | JijMatrix x -> x
    override this.ToString() =
        let jij = this.Unapply
        let sb = System.Text.StringBuilder()

        ignore <| sb.AppendLine("Ising2DContext : \n{")
        for y in 0 .. (L - 1) do
            for x in 0 .. (L - 1) do
                ignore <| sb.Append($"\t ({x}, {y}) -> [ ")
                for p in 0..3 do
                    let (a, b) = neighbour(x, y, p)
                    let f = jij.[x,y,p]
                    ignore <| sb.Append($"({a}, {b}) [{f}]")
                ignore <| sb.AppendLine($"] ")
        ignore <| sb.AppendLine("}")

        sb.ToString()

    static member ConstantFerromagnet =
        Array3D.init<float<E>> L L 4 (fun x y -> function
            | 0 -> 1.0<E> //((x + 1) % L), (y)
            | 1 -> 1.0<E> //((x + L - 1) % L), (y)
            | 2 -> 1.0<E> //(x), ((y + 1) % L)
            | 3 -> 1.0<E> //(x), ((y + L - 1) % L)
            | _ -> failwith "how?")
        |> JijMatrix

// This is the Ising2D model specification for dimension LxL
// It is bound to the interaction matrix implicitly by the size of the matrix, which makes it strictly type-unsafe, so we hide inside the module!
type private Ising2D (interactionMatrix: JijMatrix) = class
    let Jij = interactionMatrix.Unapply
    let mutable Spins = Array2D.init<bool> L L (fun _ _ -> match System.Random.Shared.Next(2) with | 0 -> false | _ -> true)
    let mutable H = 0.0<E>

    let interactionEnergy x y flip : float<E> =
        let spin' = Spins[x, y]
        let spin = if flip then not spin' else spin'

        let mutable energy = 0.0<E>
        for p in 0..3 do
            let (nx, ny) = neighbour(x, y, p)
            let neighbour_spin = Spins[nx, ny]
            let j = Jij.[x, y, p]
            energy <- if spin = neighbour_spin then energy - j else energy + j
        energy

    do
        for y in 0 .. (L - 1) do
            for x in 0 .. (L - 1) do
                H <- H + interactionEnergy x y false

    member this.Solve(T : float<K>, num_iterations : int) =
        let before = this.ToString()

        let mutable flips = 0
        let evolve x y =
            let before = interactionEnergy x y false
            let after  = interactionEnergy x y true 
            let dE = after - before
            let accept =
                if dE <= 0.0<E>
                then true
                else
                    do flips <- flips + 1
                    let beta_delta_e = dE/T
                    let probability = exp(-float(beta_delta_e))
                    let random = System.Random.Shared.NextDouble()
                    probability >= random

            if accept then
                Spins.[x, y] <- not Spins.[x, y]
                H <- H + dE

        let start = System.DateTime.Now
        for _ in 0..num_iterations do
            let (x, y) = (System.Random.Shared.Next(L), System.Random.Shared.Next(L))
            evolve x y
        let finish = System.DateTime.Now
        (before, this.ToString(), finish - start, flips)

    override _.ToString() =
        let sb = System.Text.StringBuilder()
        ignore <| sb.AppendLine("Ising2D : \n{")
        ignore <| sb.AppendLine("spins : \n\t{")
        for y in 0 .. (L - 1) do
            ignore <| sb.Append("\t\t[ ")
            for x in 0 .. (L - 1) do
                ignore <| sb.Append (if Spins[x, y] then "+ " else "0 ")
            ignore <| sb.AppendLine("]")
        ignore <| sb.AppendLine("\t}")
        ignore <| sb.AppendLine($"ham : {H}")
        ignore <| sb.AppendLine("}")

        sb.ToString()
end

let Solve (T : float<K>) (num_iterations : int) =
    JijMatrix.ConstantFerromagnet
    |> Ising2D
    |> (fun model -> model.Solve(T, num_iterations))
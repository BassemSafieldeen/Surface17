namespace Microsoft.Research.Liquid

module UserSample =
    open System
    open Util
    open Operations
    open System.Runtime.CompilerServices
    open HamiltonianGates

    //open Native             // Support for Native Interop
    //open HamiltonianGates   // Extra gates for doing Hamiltonian simulations
    //open Tests              // All the built-in tests

    /// <summary>
    /// Performs an arbitrary rotation around X. 
    /// </summary>
    /// <param name="theta">Angle to rotate by</param>
    /// <param name="qs">The head qubit of this list is operated on.</param>
    let rotX (theta:float) (qs:Qubits) =
        let gate (theta:float) =
            let nam     = "Rx" + theta.ToString("F2")
            new Gate(
                Name    = nam,
                Help    = sprintf "Rotate in X by: %f" theta,
                Mat     = (
                    let phi     = theta / 2.0
                    let c       = Math.Cos phi
                    let s       = Math.Sin phi
                    CSMat(2,[0,0,c,0.;0,1,0.,-s;1,0,0.,-s;1,1,c,0.])),
                Draw    = "\\gate{" + nam + "}"
                )
        (gate theta).Run qs

    let qfunc (qs:Qubits) =
        rotX (Math.PI/4.) qs
        for q in qs.Tail do CNOT [qs.Head;q]
        M >< qs

    [<LQD>]
    let __UserSample(n:int) =
         let stats      = Array.create 2 0
         let k          = Ket(n)
         let circ       = Circuit.Compile qfunc k.Qubits
         show "Test 1"
         circ.Dump()
         circ.RenderHT("Test 1")
         let circ       = circ.GrowGates(k)
         show "Test 2"
         circ.Dump()
         circ.RenderHT("Test 2")
         for i in 0..9999 do 
             let qs     = k.Reset(n)
             circ.Run qs
             let v      = qs.Head.Bit.v
             stats.[v] <- stats.[v] + 1
             for q in qs.Tail do
                if q.Bit <> qs.Head.Bit then
                    failwith "bad"
         show "Measured: zeros=%d and ones=%d" stats.[0] stats.[1]


    let Stabilize(qs:Qubits) =
 //     Apply the first measurement qubit
        H qs
        CNOT [qs.[0];qs.[2]]
        CNOT [qs.[0];qs.[3]]
        H qs
    //    rotX (Math.PI/2.) qs
        M qs
 //     Apply the second measurement qubit
        CNOT [qs.[1];qs.[4]]
        CNOT [qs.[7];qs.[4]]
        M qs.[4..4]
 //     Apply the third measurement qubit
        H qs.[5..5]
        CNOT [qs.[5];qs.[2]]
        CNOT [qs.[5];qs.[8]]
        CNOT [qs.[5];qs.[1]]
        CNOT [qs.[5];qs.[6]]
        H qs.[5..5]
    //    rotX (Math.PI/2.) qs.[5..5]
        M qs.[5..5]
 //     Apply the fourth measurement qubit
        CNOT [qs.[3];qs.[6]]
        CNOT [qs.[9];qs.[6]]
        CNOT [qs.[2];qs.[6]]
        CNOT [qs.[8];qs.[6]]
        M qs.[6..6]
 //     Apply the fifth measurement qubit
        CNOT [qs.[8];qs.[10]]
        CNOT [qs.[14];qs.[10]]
        CNOT [qs.[7];qs.[10]]
        CNOT [qs.[13];qs.[10]]
        M qs.[10..10]
 //     Apply the sixth measurement qubit
        H qs.[11..11]
        CNOT [qs.[11];qs.[9]]
        CNOT [qs.[11];qs.[15]]
        CNOT [qs.[11];qs.[8]]
        CNOT [qs.[11];qs.[14]]
        H qs.[11..11]
    //    rotX (Math.PI/2.) qs.[11..11]
        M qs.[11..11]
 //     Apply the seventh measurement qubit
        CNOT [qs.[9];qs.[12]]
        CNOT [qs.[15];qs.[12]]
        M qs.[12..12]
 //     Apply the eightth measurement qubit
        H qs.[16..16]
        CNOT [qs.[16];qs.[14]]
        CNOT [qs.[16];qs.[13]]
        H qs.[16..16]
    //    rotX (Math.PI/2.) qs.[16..16]
        M qs.[16..16]

(*    let test(qs:Qubits) =
        X qs
        M qs.[0..0]
        M qs.[1..1]
        M qs.[2..2]   *)

    let Stabilize2(qs:Qubits) =
        //Xs
        H qs.[0..0]
        H qs.[5..5]
        H qs.[11..11]
        H qs.[16..16]
        //top right
        CNOT [qs.[5]; qs.[2]]
        CNOT [qs.[11]; qs.[9]]
        CNOT [qs.[16]; qs.[14]]
        //top left
        CNOT [qs.[5]; qs.[1]]
        CNOT [qs.[11]; qs.[8]]
        CNOT [qs.[16]; qs.[13]]
        //bottom right
        CNOT [qs.[5]; qs.[8]]
        CNOT [qs.[11]; qs.[15]]
        CNOT [qs.[0]; qs.[3]]
        //bottom left
        CNOT [qs.[5]; qs.[7]]
        CNOT [qs.[11]; qs.[14]]
        CNOT [qs.[0]; qs.[2]]

        H qs.[0..0]
        H qs.[5..5]
        H qs.[11..11]
        H qs.[16..16]

        //Zs
        //top right
        CNOT [qs.[1]; qs.[4]]
        CNOT [qs.[3]; qs.[6]]
        CNOT [qs.[8]; qs.[10]]
        //bottom right
        CNOT [qs.[7]; qs.[5]]
        CNOT [qs.[9]; qs.[6]]
        CNOT [qs.[14]; qs.[10]]
        //top left
        CNOT [qs.[12]; qs.[9]]
        CNOT [qs.[10]; qs.[7]]
        CNOT [qs.[6]; qs.[2]]
        //bottom left
        CNOT [qs.[8]; qs.[6]]
        CNOT [qs.[13]; qs.[10]]
        CNOT [qs.[15]; qs.[12]]

        M qs.[0..0]
        M qs.[5..5]
        M qs.[11..11]
        M qs.[16..16]
        M qs.[4..4]
        M qs.[6..6]
        M qs.[10..10]
        M qs.[12..12]

    let Stabilize3(qs:Qubits) =
        Rpauli (-Math.PI/2.) Y >< qs
        //Xs
        //top right
        CZ [qs.[2]; qs.[5]]
        CZ [qs.[9]; qs.[11]]
        CZ [qs.[14]; qs.[16]]
        //top left
        CZ [qs.[1]; qs.[5]]
        CZ [qs.[8]; qs.[11]]
        CZ [qs.[13]; qs.[16]]
        //bottom right
        CZ [qs.[8]; qs.[5]]
        CZ [qs.[15]; qs.[11]]
        CZ [qs.[3]; qs.[0]]
        //bottom left
        CZ [qs.[7]; qs.[5]]
        CZ [qs.[14]; qs.[11]]
        CZ [qs.[2]; qs.[0]]

        //Zs
        //top right
        CZ [qs.[1]; qs.[4]]
        CZ [qs.[3]; qs.[6]]
        CZ [qs.[8]; qs.[10]]
        //bottom right
        CZ [qs.[7]; qs.[5]]
        CZ [qs.[9]; qs.[6]]
        CZ [qs.[14]; qs.[10]]
        //top left
        CZ [qs.[12]; qs.[9]]
        CZ [qs.[10]; qs.[7]]
        CZ [qs.[6]; qs.[2]]
        //bottom left
        CZ [qs.[8]; qs.[6]]
        CZ [qs.[13]; qs.[10]]
        CZ [qs.[15]; qs.[12]]

        Rpauli (Math.PI/2.) Y >< qs


        M qs.[0..0]
        M qs.[5..5]
        M qs.[11..11]
        M qs.[16..16]
        M qs.[4..4]
        M qs.[6..6]
        M qs.[10..10]
        M qs.[12..12]



    [<LQD>]
    let __Surface_17() =
        let surface           = Ket(17).Reset(17)
       // Rpauli (Math.PI/4.) Z surface.[8..8]; CNOT [surface.[8]; surface.[6]]; CNOT [surface.[8]; surface.[10]]; SWAP [surface.[2]; surface.[6]]; SWAP [surface.[10]; surface.[14]];  //trying state injection (could be working)
        let circ        = Circuit.Compile Stabilize3 surface
        circ.Dump()
        circ.RenderHT("Test")
(*        let circ2       = Circuit.Compile test surface
        for i in 0..999 do
            circ2.Run surface
            show "test: %d %d %d" surface.[0].Bit.v surface.[1].Bit.v surface.[2].Bit.v
            Restore [surface.[0]]
            Restore [surface.[1]]
            Restore [surface.[2]] *)

        for i in 0..999 do
            circ.Run surface
            show "Syndrome measurements: %d %d %d %d %d %d %d %d" surface.[0].Bit.v surface.[4].Bit.v surface.[5].Bit.v surface.[6].Bit.v surface.[10].Bit.v surface.[11].Bit.v surface.[12].Bit.v surface.[16].Bit.v
            Reset Zero [surface.[0]]; Reset Zero [surface.[4]]; Reset Zero [surface.[5]]; Reset Zero [surface.[6]]; Reset Zero [surface.[10]]; Reset Zero [surface.[11]]; Reset Zero [surface.[12]]; Reset Zero [surface.[16]];
            if i=10 then
                show "Logical H"
                H surface.[1..1]; H surface.[2..2]; H surface.[3..3]; H surface.[7..7]; H surface.[8..8]; H surface.[9..9]; H surface.[13..13]; H surface.[14..14]; H surface.[15..15];
                
                //show "Logical Z"
                //Z surface.[2..2]; Z surface.[8..8]; Z surface.[14..14];

                //show "Logical X"
                //X surface.[7..7]; X surface.[8..8]; X surface.[9..9];

module Main =
    open App

    /// <summary>
    /// The main entry point for Liquid.
    /// </summary>
    [<EntryPoint>]
    let Main _ =
        RunLiquid ()

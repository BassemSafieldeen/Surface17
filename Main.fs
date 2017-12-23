namespace Microsoft.Research.Liquid

module UserSample =
    open System
    open Util
    open Operations
    open System.Runtime.CompilerServices
    open HamiltonianGates
    open System.Text
    //open Native             // Support for Native Interop
    //open HamiltonianGates   // Extra gates for doing Hamiltonian simulations
    //open Tests              // All the built-in tests

    /// <summary>
    /// Performs an arbitrary rotation around X. 
    /// </summary>
    /// <param name="theta">Angle to rotate by</param>
    /// <param name="qs">The head qubit of this list is operated on.</param>

    //================= following are two examples : QECC() and NoiseAmpp() =============
    let myEPR (qs:Qubits) = 
        H qs; CNOT qs;  //X qs; // M qs
    let myfun (qs:Qubits) =
        H qs; M >< qs       
    [<LQD>]
    let QECC() =  // A smiple example showing how Stabilizer and its methods work. This code is derived from Liquid/Sample/QECC.fsx 
        let k               = Ket(2)
        let qs               = k.Reset(2)
        let circc            = Circuit.Compile myEPR qs  // create a circuit first
        let stab            = Stabilizer(circc,k)        // generate a Stabilizer object
        
        stab.Run()
        //let _,b0            = stab.[0]    // If there's "M qs[0]" in circc, stab.[0] would show the measurement result of qs[0]. If not, then this line would cause run time error.
        //let _,b1            = stab.[1]
        //show "EPR in stabilizer: [%d%d] " b0.v b1.v
        show ""
        // Show the final state in the form of stabilizer tableau. 
        show "=== Final State: "  
        stab.ShowState showInd 0
        stab.Gaussian()
        show "=== After Gaussian: "
        stab.ShowState showInd 0
        show ""
        (*  // This comment shows how to intepret the screen print. Tkae the final state of an EPR pair as an example :
        0:0000.0/=== Final State:
        0:0000.0/
        0:0000.0/+Z.
        0:0000.0/+.X
        0:0000.0/---
        0:0000.0/+XX    <----------- The final state is stablized by X1X2. That is, the eigenstate of X1X2 with eigenvalue +1.
        0:0000.0/+ZZ    <-----------                              by Z1Z2.                            Z1Z2                 +1.

        The upper/lower part is destabilizer/stabilizer generators. see "An Introduction to Stabilizer Circuit Simulation UMD Department of ..."
        
        In this EPR example, Z1 and X2 are destabilizer generators. X1X2, Z1Z2 are stabilizer generators.
        *)
    [<LQD>]
    let NoiseAmpp() =  // A example code showing how to obtain the state without doing state tomography. This code is derived from Liquid/Sample/NoiseAmp.fsx 
                       // This may be helpful in understanding NoiseAmp.fsx.
        // Output dump routine
        let sb              = StringBuilder()
        let app (x:string)  = sb.Append x |> ignore
        let dump (m:bool) (iter:int) (v:CVec) =      //show every component of ket v
            if iter = 0 then show "Iter,qs=00,qs=01,qs=10,qs=11"

            sb.Length      <- 0
            sprintf "%4d" iter |> app
            for i in 0UL..v.Length-1UL do
                app ","
                if m = true then
                    sprintf "%7.3f" v.[i].r |> app  // get the real part of v.[i]
                    sprintf "+%7.3f i" v.[i].i |> app  // get the real part of v.[i]
                if m = false then
                    sprintf "%7.5f" v.[i].MCC |> app  // get the real part of v.[i]
            show "%O" sb

        // 2 Qubit tests
        let ket     = Ket(2)
        let qs          = ket.Reset(2)    
        let circ    = Circuit.Compile (fun (qs:Qubits) -> Rpauli (Math.PI/8.) X  qs) ket.Qubits
        
        //Get a handle to the state vector for output
        let v           = ket.Single()  // fully realized state vector (2^n in size)
        dump true 0 v
        for iter in 1..30 do
            circ.Run qs 
            dump true iter v  //show the state ket = a|00> + b|01> + c|10> + d|11>
        
        let qs          = ket.Reset(2) 
        let v           = ket.Single()
        dump false 0 v
        for iter in 1..30 do
            circ.Run qs 
            dump false iter v  //show magnitude square of ket (probability of being 0 and 1 )

      //================= above are two examples : QECC() and NoiseAmpp() =============

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
        M qs.[2..2]   *)  // Comment

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

        
    let Stabilize4(qs:Qubits) =
        H >< qs
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
        
        H qs.[0..0];        H qs.[1..1];        H qs.[2..2];        H qs.[3..3]  // 4
        H qs.[5..5]       // 6
        H qs.[7..7];        H qs.[8..8];        H qs.[9..9]        // 10
        H qs.[11..11]        // 12
        H qs.[13..13];        H qs.[14..14];        H qs.[15..15];        H qs.[16..16]
        //Zs
        //top right
        CZ [qs.[1]; qs.[4]]
        CZ [qs.[3]; qs.[6]]
        CZ [qs.[8]; qs.[10]]
        //bottom right
        CZ [qs.[7]; qs.[4]]
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

        H qs.[4..4];         H qs.[6..6];        H qs.[10..10];        H qs.[12..12]


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
        let stats                  = Array.create 2 0      
        for j in 0..100 do  //this loop is here because I am testing simple tomography on the decoded logical state
            let surface           = Ket(17).Reset(17)
            Rpauli (Math.PI/4.) Y surface.[8..8]; CNOT [surface.[8]; surface.[6]]; CNOT [surface.[8]; surface.[10]]; SWAP [surface.[2]; surface.[6]]; SWAP [surface.[10]; surface.[14]];  //trying state injection (could be working)

      (*      let stat              = Array.create 2 0
            for i in 0..50 do
                Rpauli (Math.PI/4.) Y surface.[8..8]; //CNOT [surface.[8]; surface.[6]]; CNOT [surface.[8]; surface.[10]]; SWAP [surface.[2]; surface.[6]]; SWAP [surface.[10]; surface.[14]];  //trying state injection (could be working)
                M surface.[8..8]
                stat.[0 + surface.[8].Bit.v] <- stat.[0 + surface.[8].Bit.v] + 1 
                show "stats: Zeros=%d Ones=%d" stat.[0] stat.[1]
                Reset Zero [surface.[8]];   *)
            let circ        = Circuit.Compile Stabilize4 surface
            circ.Dump()
            circ.RenderHT("Test")
            
    (*        let circ2       = Circuit.Compile test surface
            for i in 0..999 do
                circ2.Run surface
                show "test: %d %d %d" surface.[0].Bit.v surface.[1].Bit.v surface.[2].Bit.v
                Restore [surface.[0]]
                Restore [surface.[1]]
                Restore [surface.[2]] *)

            for i in 0..5 do
                circ.Run surface
                if i <> 5 then   //this condition is here because at the last step we need to do decoding
                    show "Syndrome measurements: %d %d %d %d %d %d %d %d" surface.[0].Bit.v surface.[4].Bit.v surface.[5].Bit.v surface.[6].Bit.v surface.[10].Bit.v surface.[11].Bit.v surface.[12].Bit.v surface.[16].Bit.v
                    Reset Zero [surface.[0]]; Reset Zero [surface.[4]]; Reset Zero [surface.[5]]; Reset Zero [surface.[6]]; Reset Zero [surface.[10]]; Reset Zero [surface.[11]]; Reset Zero [surface.[12]]; Reset Zero [surface.[16]];
                if i=2 then
                    show "__"

                    //show "Logical H"
                    //H surface.[1..1]; H surface.[2..2]; H surface.[3..3]; H surface.[7..7]; H surface.[8..8]; H surface.[9..9]; H surface.[13..13]; H surface.[14..14]; H surface.[15..15];
                
                    //show "Logical Z"
                    //Z surface.[2..2]; Z surface.[8..8]; Z surface.[14..14];

                    //show "Logical X"
                    //X surface.[7..7]; X surface.[8..8]; X surface.[9..9];

            //Decoding the logical State
            show "Decoding the Logical State"
            //show "MA3 = %d and MA4 = %d" surface.[6].Bit.v surface.[10].Bit.v
            M surface.[7..7]; M surface.[9..9];
            show "MD3 = %d and MD5 = %d" surface.[7].Bit.v surface.[9].Bit.v
            if ((surface.[7].Bit.v + surface.[9].Bit.v)%2) = 0 then
                show "No X needed"
            else 
                show "X needed"
                X surface.[8..8]

            CNOT [surface.[8]; surface.[2]]; CNOT [surface.[8]; surface.[14]];// SWAP [surface.[2]; surface.[6]]; SWAP [surface.[10]; surface.[14]]; 
            M surface.[1..1]; M surface.[2..2]; M surface.[3..3]; M surface.[13..13]; M surface.[14..14]; M surface.[15..15];
            //collapse all the other qubits. If there's entanglement between surf[8] then this would differ from just measure surf[8]
            
            //Rpauli (-Math.PI/2.) Y surface.[8..8];    //state tomo X
            //Rpauli (Math.PI/2.) X surface.[8..8];     //state tomo Y
            
            M surface.[8..8];
            stats.[0 + surface.[8].Bit.v] <- stats.[0 + surface.[8].Bit.v] + 1; 
            show "stats: Zeros = %d   and Ones = %d" stats.[0] stats.[1]

module Main =
    open App

    /// <summary>
    /// The main entry point for Liquid.
    /// </summary>
    [<EntryPoint>]
    let Main _ =
        RunLiquid ()

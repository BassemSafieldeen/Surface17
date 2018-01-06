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
    let QECC() =  // A simple example showing how Stabilizer and its methods work. This code is derived from Liquid/Sample/QECC.fsx 
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
        let v           = ket.Single() // fully realized state vector (2^n in size)
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

    let sb              = StringBuilder()
    let app (x:string)  = sb.Append x |> ignore
    let dump (m:bool) (iter:int) (v:CVec) =      //show every component of ket v
        //if iter = 0 then show "Iter,qs=00,qs=01,qs=10,qs=11"

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

    [<LQD>]
    let __Surface_17() =
        let stats                  = Array.create 2 0      
        let ket               = Ket(17)
        let surface           = ket.Reset(17)
        //State Injection
        Rpauli (Math.PI/8.) Y surface.[8..8]; 
        //let v = ket.Single()    // Uncomment this line and the following line to see the state after injection
        //dump false 0 v
        CNOT [surface.[8]; surface.[6]]; CNOT [surface.[8]; surface.[10]]; SWAP [surface.[2]; surface.[6]]; SWAP [surface.[10]; surface.[14]];  //State Injection
        let circ        = Circuit.Compile Stabilize ket.Qubits
        circ.Dump()
        circ.RenderHT("Test")


        for i in 0..5 do
            circ.Run surface
            if i <> 5 then   //this condition is here because at the last step we need to do decoding
                show "Syndrome measurements: %d %d %d %d %d %d %d %d" surface.[0].Bit.v surface.[4].Bit.v surface.[5].Bit.v surface.[6].Bit.v surface.[10].Bit.v surface.[11].Bit.v surface.[12].Bit.v surface.[16].Bit.v
                Reset Zero [surface.[0]]; Reset Zero [surface.[4]]; Reset Zero [surface.[5]]; Reset Zero [surface.[6]]; Reset Zero [surface.[10]]; Reset Zero [surface.[11]]; Reset Zero [surface.[12]]; Reset Zero [surface.[16]];
            //if i=2 then
                //show "__"

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


        //State Tomography
        show "Doing Tomography (this destroys the surface, i.e. you cannot do tomography and carry on with other operations)."
        for i in 0..16 do
            if i <> 8 then
                Reset Zero [surface.[i]];
        let v = ket.Single()
        dump false 0 v


module Main =
    open App

    /// <summary>
    /// The main entry point for Liquid.
    /// </summary>
    [<EntryPoint>]
    let Main _ =
        RunLiquid ()

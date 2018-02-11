﻿namespace Microsoft.Research.Liquid

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
        CZ [qs.[1]; qs.[11]]
        CZ [qs.[5]; qs.[14]]
        CZ [qs.[7]; qs.[16]]
        //top left
        CZ [qs.[0]; qs.[11]]
        CZ [qs.[4]; qs.[14]]
        CZ [qs.[6]; qs.[16]]
        //bottom right
        CZ [qs.[4]; qs.[11]]
        CZ [qs.[8]; qs.[14]]
        CZ [qs.[2]; qs.[9]]
        //bottom left
        CZ [qs.[3]; qs.[11]]
        CZ [qs.[7]; qs.[14]]
        CZ [qs.[1]; qs.[9]]
        
        H qs.[9..9];        H qs.[0..0];        H qs.[1..1];        H qs.[2..2]  // 4
        H qs.[11..11]       // 6
        H qs.[3..3];        H qs.[4..4];        H qs.[5..5]        // 10
        H qs.[14..14]        // 12
        H qs.[6..6];        H qs.[7..7];        H qs.[8..8];        H qs.[16..16]
        //Zs
        //top right
        CZ [qs.[0]; qs.[10]]
        CZ [qs.[2]; qs.[12]]
        CZ [qs.[4]; qs.[13]]
        //bottom right
        CZ [qs.[3]; qs.[10]]
        CZ [qs.[5]; qs.[12]]
        CZ [qs.[7]; qs.[13]]
        //top left
        CZ [qs.[15]; qs.[5]]
        CZ [qs.[13]; qs.[3]]
        CZ [qs.[12]; qs.[1]]
        //bottom left
        CZ [qs.[4]; qs.[12]]
        CZ [qs.[6]; qs.[13]]
        CZ [qs.[8]; qs.[15]]

        H qs.[10..10];         H qs.[12..12];        H qs.[13..13];        H qs.[15..15]


        M qs.[9..9]
        M qs.[11..11]
        M qs.[14..14]
        M qs.[16..16]
        M qs.[10..10]
        M qs.[12..12]
        M qs.[13..13]
        M qs.[15..15]

    let sb              = StringBuilder()
    let app (x:string)  = sb.Append x |> ignore
    let dump (m:bool) (iter:int) (v:CVec) =      //show every component of ket v
        //if iter = 0 then show "Iter,qs=00,qs=01,qs=10,qs=11"

        sb.Length      <- 0
        //sprintf "%4d" iter |> app
        for i in 0UL..v.Length-1UL do
            //app ", "
            if m = true then
                if v.[i].r <> 0.0 || v.[i].i <> 0.0 then 
                    sprintf "  C_%i = " i |> app
                    if v.[i].r <> 0.0 then sprintf "%4.3f" v.[i].r |> app  // get the real part of v.[i]
                    if v.[i].i <> 0.0 then
                        if v.[i].i > 0.0 then sprintf "+%4.3f i" v.[i].i |> app  // get the real part of v.[i]
                        else sprintf "%4.3f i" v.[i].i |> app  
                    sprintf "  ,  " |> app
            if m = false then
                if v.[i].MCC <> 0.0 then sprintf ", i=%4d" i |> app; sprintf " MCC=%7.5f" v.[i].MCC |> app  // get the real part of v.[i]
        show "%O" sb

    [<LQD>]
    let __Surface_17() =
        let stats                  = Array.create 2 0 
        let mutable prev_syndrome =  [|0; 0; 0; 0; 0; 0; 0; 0|];
        let mutable new_syndrome = [|0; 0; 0; 0; 0; 0; 0; 0|];
        let mutable changes = [|0; 0; 0; 0; 0; 0; 0; 0|];
        let ket               = Ket(17);
        let surface           = ket.Reset(17);
        let mutable flag = 0;


        let theta           = Math.PI*0.89        //specify the single qubit state to be injected
        let phi             = Math.PI*1.56        //specify the single qubit state to be injected
        let probDamp        = 0.//2.e-2             //amplitude damping noise
        let probPolar       = 1.e-2         //depolaried noise
        let N_cycle         = 50

        //Prepare arbitrary single qubit state
        Rpauli (-theta) Y surface.[4..4]; 
        Rpauli (phi) Z surface.[4..4];
        let v = ket.Single()    // Uncomment this line and the following line to see the state after injection
        show "  Initial state of data qubit 4 = "
        dump true 0 v
        Console.ReadLine() |> ignore
        //State Injection
        CNOT [surface.[4]; surface.[12]]; CNOT [surface.[4]; surface.[13]]; SWAP [surface.[1]; surface.[12]]; SWAP [surface.[13]; surface.[7]];  //State Injection
        let circ        = Circuit.Compile Stabilize ket.Qubits
        circ.Dump()
        circ.RenderHT("Test")

        // Create noise model
        // Probabilities for our two types of noise
        let circN    = Circuit.Compile (fun (qs:Qubits) -> I qs.[0..0];I qs.[1..1]; I qs.[2..2];I qs.[3..3];I qs.[4..4];I qs.[5..5];I qs.[6..6];I qs.[7..7];I qs.[8..8]; ) ket.Qubits //noise channel
        

        let mkM (p:float) (g:string) (mx:int) = {Noise.DefaultNoise p with gate=g;maxQs=mx}
        let models      = [  mkM probPolar "I" 1 ]

        let noise           = Noise(circN,ket,models)
        
        ket.TraceRun       <- 0         // 1=log 2=console
        noise.LogGates     <- false     // Show each gate execute?
        noise.TraceWrap    <- false
        noise.TraceNoise   <- false
        noise.DampProb(0)  <- probDamp ; noise.DampProb(1)  <- probDamp; noise.DampProb(2)  <- probDamp;// apply damping error on qubit 1
        noise.DampProb(3)  <- probDamp ; noise.DampProb(4)  <- probDamp; noise.DampProb(5)  <- probDamp;
        noise.DampProb(6)  <- probDamp ; noise.DampProb(7)  <- probDamp; noise.DampProb(8)  <- probDamp;             
        // End noise model


        for i in 0..N_cycle do
            if i <>0 then noise.Run ket
            circ.Run ket.Qubits
            noise.Dump(showInd,0,true)
            if i <> N_cycle then   //this condition is here because at the last step we need to do decoding
                show " At round %i, Syndrome measurements: %d %d %d %d %d %d %d %d" i surface.[9].Bit.v surface.[10].Bit.v surface.[11].Bit.v surface.[12].Bit.v surface.[13].Bit.v surface.[14].Bit.v surface.[15].Bit.v surface.[16].Bit.v
                //Implementing Decoder
                new_syndrome <- [|surface.[9].Bit.v; surface.[10].Bit.v; surface.[11].Bit.v; surface.[12].Bit.v; surface.[13].Bit.v; surface.[14].Bit.v; surface.[15].Bit.v; surface.[16].Bit.v |];
                if i = 0 then prev_syndrome <- new_syndrome;
                changes <- [| for k in 0..7 -> prev_syndrome.[k] ^^^ new_syndrome.[k]|]; //^^^ is an xor
                if (changes.[0] = 1) && (changes.[2] = 0) then show "Z error on data qubit 2 occurred. Applying fix."; Z ket.Qubits.[2..2]; flag <- 1;
                if (changes.[0] = 1) && (changes.[2] = 1) then show "Z error on data qubit 1 occurred. Applying fix."; Z ket.Qubits.[1..1]; flag <- 1;
                if (changes.[1] = 1) && (changes.[4] = 0) then show "X error on data qubit 0 occurred. Applying fix."; X ket.Qubits.[0..0]; flag <- 1;
                if (changes.[1] = 1) && (changes.[4] = 1) then show "X error on data qubit 3 occurred. Applying fix."; X ket.Qubits.[3..3]; flag <- 1;
                if (changes.[2] = 1) && (changes.[0] = 0) && (changes.[5] = 0) then show "Z error on data qubits 3 or 0 occurred. Applying fix."; Z ket.Qubits.[0..0]; flag <- 1;
                //if (changes.[2] = 1) && (changes.[0] = 1) && (changes.[5] = 0) then show "Z error on data qubit 1 occurred. Applying fix."; Z ket.Qubits.[1..1]; flag <- 1; redundant
                if (changes.[2] = 1) && (changes.[0] = 0) && (changes.[5] = 1) then show "Z error on data qubit 4 occurred. Applying fix."; Z ket.Qubits.[4..4]; flag <- 1;
                if (changes.[3] = 1) && (changes.[4] = 0) && (changes.[6] = 0) then show "X error on data qubits 1 or 2 occurred. Applying fix."; X ket.Qubits.[1..1]; flag <- 1;
                if (changes.[3] = 1) && (changes.[4] = 1) && (changes.[6] = 0) then show "X error on data qubit 4 occurred. Applying fix."; X ket.Qubits.[4..4]; flag <- 1;
                //if (changes.[3] = 1) && (changes.[4] = 0) && (changes.[6] = 1) then show "X error on data qubit 5 occurred. Applying fix."; X ket.Qubits.[5..5]; flag <- 1; redundant
                if (changes.[4] = 1) && (changes.[1] = 0) && (changes.[3] = 0) then show "X error on data qubits 6 or 7 occurred. Applying fix."; X ket.Qubits.[6..6]; flag <- 1;
                //if (changes.[4] = 1) && (changes.[1] = 1) && (changes.[3] = 0) then show "X error on data qubit 3 occurred"; redundant
                //if (changes.[4] = 1) && (changes.[1] = 0) && (changes.[3] = 1) then show "X error on data qubit 4 occurred"; redundant
                if (changes.[5] = 1) && (changes.[2] = 0) && (changes.[7] = 0) then show "Z error on data qubits 5 or 8 occurred. Applying fix."; Z ket.Qubits.[8..8]; flag <- 1;
                //if (changes.[5] = 1) && (changes.[2] = 1) && (changes.[7] = 0) then show "Z error on data qubit 4 occurred"; redundant
                if (changes.[5] = 1) && (changes.[2] = 0) && (changes.[7] = 1) then show "Z error on data qubit 7 occurred. Applying fix."; Z ket.Qubits.[7..7]; flag <- 1;
                if (changes.[6] = 1) && (changes.[3] = 0) then show "X error on data qubit 8 occurred. Applying fix."; X ket.Qubits.[8..8]; flag <- 1;
                if (changes.[6] = 1) && (changes.[3] = 1) then show "X error on data qubit 5 occurred. Applying fix."; X ket.Qubits.[5..5]; flag <- 1;
                if (changes.[7] = 1) && (changes.[5] = 0) then show "Z error on data qubit 6 occurred. Applying fix."; Z ket.Qubits.[6..6]; flag <- 1;
                //if (changes.[7] = 1) && (changes.[5] = 1) then show "Z error on data qubit 7 occurred"; redundant
                if flag = 0 then prev_syndrome <- new_syndrome;
                if flag = 1 then Console.ReadLine() |> ignore
                flag <- 0;
                
                Reset Zero [surface.[9]]; Reset Zero [surface.[10]]; Reset Zero [surface.[11]]; Reset Zero [surface.[12]]; Reset Zero [surface.[13]]; Reset Zero [surface.[14]]; Reset Zero [surface.[15]]; Reset Zero [surface.[16]];
            //if i=2 then
                //show "__"

                //show "Logical H"
                //H surface.[0..0]; H surface.[1..1]; H surface.[2..2]; H surface.[3..3]; H surface.[4..4]; H surface.[5..5]; H surface.[6..6]; H surface.[7..7]; H surface.[8..8];
                
                //show "Logical Z"
                //Z surface.[0..0]; Z surface.[4..4]; Z surface.[8..8];

                //show "Logical X"
                //X surface.[2..2]; X surface.[4..4]; X surface.[6..6];

        //Decoding the logical State
        show "Decoding the Logical State"
        //show "MA3 = %d and MA4 = %d" surface.[12].Bit.v        ` surface.[13].Bit.v
        M surface.[3..3]; M surface.[5..5];
        show "MD3 = %d and MD5 = %d" surface.[3].Bit.v surface.[5].Bit.v
        if ((surface.[3].Bit.v + surface.[5].Bit.v)%2) = 0 then
            show "No X needed"
        else 
            show "X needed"
            X surface.[4..4]
        CNOT [surface.[4]; surface.[1]]; CNOT [surface.[4]; surface.[7]];// SWAP [surface.[1]; surface.[12]]; SWAP [surface.[7]; surface.[13]]; 
        M surface.[0..0]; M surface.[1..1]; M surface.[2..2]; M surface.[6..6]; M surface.[7..7]; M surface.[8..8];
        //collapse all the other qubits. If there's entanglement between surf[8] then this would differ from just measure surf[8] 


        //State Tomography
        show "Doing Tomography (this destroys the surface and leaves only the central data qubit, i.e. you cannot do tomography and carry on with other operations)."
        for i in 0..16 do
            if i <> 4 then
                Reset Zero [surface.[i]];
        let v = ket.Single()
        show "  Final state of data qubit 4 = "
        dump true 0 v


module Main =
    open App

    /// <summary>
    /// The main entry point for Liquid.
    /// </summary>
    [<EntryPoint>]
    let Main _ =
        RunLiquid ()

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

    let StabilizH(qs:Qubits) =
        H >< qs
        //Xs
        //top right
        CZ [qs.[0]; qs.[10]]
        CZ [qs.[2]; qs.[12]]
        CZ [qs.[4]; qs.[13]]
        //top left
        CZ [qs.[1]; qs.[12]]
        CZ [qs.[3]; qs.[13]]
        CZ [qs.[5]; qs.[15]]
        //bottom right
        CZ [qs.[3]; qs.[10]]
        CZ [qs.[5]; qs.[12]]
        CZ [qs.[7]; qs.[13]]
        //bottom left
        CZ [qs.[4]; qs.[12]]
        CZ [qs.[6]; qs.[13]]
        CZ [qs.[8]; qs.[15]]
        
        H qs.[10..10];        H qs.[0..0];        H qs.[1..1];        H qs.[2..2]  // 4
        H qs.[12..12]       // 6
        H qs.[3..3];        H qs.[4..4];        H qs.[5..5]        // 10
        H qs.[13..13]        // 12
        H qs.[6..6];        H qs.[7..7];        H qs.[8..8];        H qs.[15..15]
        //Zs
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

        H qs.[9..9];         H qs.[11..11];        H qs.[14..14];        H qs.[16..16]


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
        let mutable prev_syndrome =  [|0; 0; 0; 0; 0; 0; 0; 0|];
        let mutable new_syndrome = [|0; 0; 0; 0; 0; 0; 0; 0|];
        let mutable changes = [|0; 0; 0; 0; 0; 0; 0; 0|];
        let ket               = Ket(17);
        let surface           = ket.Reset(17);
        let mutable flag = 0;
        let mutable flagH = 0;
        let mutable sZ = "Z";
        let mutable sX = "X";

        let theta           = Math.PI*0.7;        //specify the single qubit state to be injected
        let phi             = Math.PI*0.;       //specify the single qubit state to be injected
        let probDamp        = 0.//1.e-2;             //amplitude damping noise
        let probPolar       = 3.e-2;         //depolaried noise
        let N_cycle         = 25;
        let apply_logical_X = 60;
        let apply_logical_Z = 60;
        let apply_logical_H = 10;

        //Prepare arbitrary single qubit state
        Rpauli (-theta) Y surface.[4..4]; 
        Rpauli (phi) Z surface.[4..4];
        let v = ket.Single();    // Uncomment this line and the following line to see the state after injection
        show "  Initial state of data qubit 4 = "
        dump true 0 v;
        Console.ReadLine() |> ignore;
        //State Injection
        CNOT [surface.[4]; surface.[12]]; CNOT [surface.[4]; surface.[13]]; SWAP [surface.[1]; surface.[12]]; SWAP [surface.[13]; surface.[7]];  //State Injection
        let circ        = Circuit.Compile Stabilize ket.Qubits;
        let circH       = Circuit.Compile StabilizH ket.Qubits;
        circ.Dump();
        circ.RenderHT("Test");

        // Create noise model
        // Probabilities for our two types of noise
        let circN    = Circuit.Compile (fun (qs:Qubits) -> I qs.[0..0];I qs.[1..1]; I qs.[2..2];I qs.[3..3];I qs.[4..4];I qs.[5..5];I qs.[6..6];I qs.[7..7];I qs.[8..8]; ) ket.Qubits //noise channel
        

        let mkM (p:float) (g:string) (mx:int) = {Noise.DefaultNoise p with gate=g;maxQs=mx}
        let models      = [  mkM probPolar "I" 1 ]

        let noise           = Noise(circN,ket,models);
        
        ket.TraceRun       <- 0;         // 1=log 2=console
        noise.LogGates     <- true;     // Show each gate execute?
        noise.TraceWrap    <- false;
        noise.TraceNoise   <- false;
        noise.DampProb(0)  <- probDamp ; noise.DampProb(1)  <- probDamp; noise.DampProb(2)  <- probDamp;// apply damping error on qubit 1
        noise.DampProb(3)  <- probDamp ; noise.DampProb(4)  <- probDamp; noise.DampProb(5)  <- probDamp;
        noise.DampProb(6)  <- probDamp ; noise.DampProb(7)  <- probDamp; noise.DampProb(8)  <- probDamp;             
        // End noise model


        for i in 0..N_cycle do
            if i <> 0 then noise.Run ket;
            if flagH = 0 then circ.Run ket.Qubits else circH.Run ket.Qubits;
            noise.Dump(showInd,0,true);
            
            show " At round %i, Syndrome measurements: %d %d %d %d %d %d %d %d" i surface.[9].Bit.v surface.[10].Bit.v surface.[11].Bit.v surface.[12].Bit.v surface.[13].Bit.v surface.[14].Bit.v surface.[15].Bit.v surface.[16].Bit.v
            //Implementing Decoder
            new_syndrome <- [|surface.[9].Bit.v; surface.[10].Bit.v; surface.[11].Bit.v; surface.[12].Bit.v; surface.[13].Bit.v; surface.[14].Bit.v; surface.[15].Bit.v; surface.[16].Bit.v |];
            if i = 0 then prev_syndrome <- new_syndrome;
            changes <- [| for k in 0..7 -> prev_syndrome.[k] ^^^ new_syndrome.[k]|]; //^^^ is an xor
            if flagH = 1 then H >< surface;          
            if (changes.[0] = 1) && (changes.[2] = 0) then show "%s error on data qubit 2 occurred. Applying fix." sZ; Z ket.Qubits.[2..2]; flag <- 1;
            if (changes.[0] = 1) && (changes.[2] = 1) then show "%s error on data qubit 1 occurred. Applying fix." sZ; Z ket.Qubits.[1..1]; flag <- 1;
            if (changes.[1] = 1) && (changes.[4] = 0) then show "%s error on data qubit 0 occurred. Applying fix." sX; X ket.Qubits.[0..0]; flag <- 1;
            if (changes.[1] = 1) && (changes.[4] = 1) then show "%s error on data qubit 3 occurred. Applying fix." sX; X ket.Qubits.[3..3]; flag <- 1;
            if (changes.[2] = 1) && (changes.[0] = 0) && (changes.[5] = 0) then show "%s error on data qubits 3 or 0 occurred. Applying fix." sZ; Z ket.Qubits.[0..0]; flag <- 1;
            if (changes.[2] = 1) && (changes.[0] = 0) && (changes.[5] = 1) then show "%s error on data qubit 4 occurred. Applying fix." sZ; Z ket.Qubits.[4..4]; flag <- 1;
            if (changes.[3] = 1) && (changes.[4] = 0) && (changes.[6] = 0) then show "%s error on data qubits 1 or 2 occurred. Applying fix." sX; X ket.Qubits.[1..1]; flag <- 1;
            if (changes.[3] = 1) && (changes.[4] = 1) && (changes.[6] = 0) then show "%s error on data qubit 4 occurred. Applying fix." sX; X ket.Qubits.[4..4]; flag <- 1;
            if (changes.[4] = 1) && (changes.[1] = 0) && (changes.[3] = 0) then show "%s error on data qubits 6 or 7 occurred. Applying fix." sX; X ket.Qubits.[6..6]; flag <- 1;
            if (changes.[5] = 1) && (changes.[2] = 0) && (changes.[7] = 0) then show "%s error on data qubits 5 or 8 occurred. Applying fix." sZ; Z ket.Qubits.[8..8]; flag <- 1;
            if (changes.[5] = 1) && (changes.[2] = 0) && (changes.[7] = 1) then show "%s error on data qubit 7 occurred. Applying fix." sZ; Z ket.Qubits.[7..7]; flag <- 1;
            if (changes.[6] = 1) && (changes.[3] = 0) then show "%s error on data qubit 8 occurred. Applying fix." sX; X ket.Qubits.[8..8]; flag <- 1;
            if (changes.[6] = 1) && (changes.[3] = 1) then show "%s error on data qubit 5 occurred. Applying fix." sX; X ket.Qubits.[5..5]; flag <- 1;
            if (changes.[7] = 1) && (changes.[5] = 0) then show "%s error on data qubit 6 occurred. Applying fix." sZ; Z ket.Qubits.[6..6]; flag <- 1;
            if flagH = 1 then H >< surface;
            if flag = 0 then prev_syndrome <- new_syndrome;
            if flag = 1 then Console.ReadLine() |> ignore;
            flag <- 0;
            if i <> N_cycle then   //this condition is here because at the last step we need to do decoding    
                Reset Zero [surface.[9]]; Reset Zero [surface.[10]]; Reset Zero [surface.[11]]; Reset Zero [surface.[12]]; Reset Zero [surface.[13]]; Reset Zero [surface.[14]]; Reset Zero [surface.[15]]; Reset Zero [surface.[16]];
            
            if i = apply_logical_H  then
                show "Apply logical H" ; Console.ReadLine() |> ignore;
                H surface.[0..0]; H surface.[1..1]; H surface.[2..2]; H surface.[3..3]; H surface.[4..4]; H surface.[5..5]; H surface.[6..6]; H surface.[7..7]; H surface.[8..8];
                if flagH = 0 then flagH <- 1; sZ <- "X"; sX <- "Z"; else flagH <- 0;sZ <- "Z"; sX <- "X";
            if i = apply_logical_Z then     
                show "Apply logical Z" ; Console.ReadLine() |> ignore;
                if flagH =0 then Z surface.[3..3]; Z surface.[4..4]; Z surface.[5..5]; else Z surface.[1..1]; Z surface.[4..4]; Z surface.[7..7];
            if i = apply_logical_X then
                show "Apply logical X" ;Console.ReadLine() |> ignore;
                if flagH = 0 then X surface.[1..1]; X surface.[4..4]; X surface.[7..7]; else X surface.[3..3]; X surface.[4..4]; X surface.[5..5];

        //Decoding the logical State
        show "Decoding the Logical State";
        if flagH = 0 then
            M surface.[3..3]; M surface.[5..5];        show "MD3 = %d and MD5 = %d" surface.[3].Bit.v surface.[5].Bit.v;
            if ((surface.[3].Bit.v + surface.[5].Bit.v)%2) = 0 then     show "No X needed";   
            else  show "X needed";  X surface.[4..4];
            CNOT [surface.[4]; surface.[1]]; CNOT [surface.[4]; surface.[7]];// SWAP [surface.[1]; surface.[12]]; SWAP [surface.[7]; surface.[13]]; 
            M surface.[0..0]; M surface.[1..1]; M surface.[2..2]; M surface.[6..6]; M surface.[7..7]; M surface.[8..8];
        //collapse all the other qubits. If there's entanglement between surf[8] then this would differ from just measure surf[8] 
        else 
            M surface.[1..1]; M surface.[7..7];        show "MD1 = %d and MD7 = %d" surface.[1].Bit.v surface.[7].Bit.v;
            if ((surface.[1].Bit.v + surface.[7].Bit.v)%2) = 0 then     show "No X needed";   
            else  show "X needed";  X surface.[4..4];
            CNOT [surface.[4]; surface.[3]]; CNOT [surface.[4]; surface.[5]];// SWAP [surface.[3]; surface.[12]]; SWAP [surface.[5]; surface.[13]]; 
            M surface.[0..0]; M surface.[2..2]; M surface.[3..3]; M surface.[5..5]; M surface.[6..6]; M surface.[8..8];

        


        //State Tomography
        show "Doing Tomography (this destroys the surface and leaves only the central data qubit, i.e. you cannot do tomography and carry on with other operations)."
        for i in 0..16 do
            if i <> 4 then
                Reset Zero [surface.[i]];
        let v = ket.Single()
        show "  Final state of data qubit 4 = ";
        dump true 0 v;


module Main =
    open App

    /// <summary>
    /// The main entry point for Liquid.
    /// </summary>
    [<EntryPoint>]
    let Main _ =
        RunLiquid ()

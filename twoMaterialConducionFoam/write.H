    if (runTime.outputTime())
    {
           volVectorField gradT = fvc::grad(T);

            runTime.write();
    }

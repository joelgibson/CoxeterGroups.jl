// The purpose of this script is to generate test-cases for CoxGrpMin_ReflTable.jl.
//
// To spot-check the minimal root implementation, we will copy some tables generated by Magma into the tests.
// These tables need to be modified a little for our Julia implementation: Magma uses refl[s, α_s] = -s wheras
// we just use zero, and Magma uses some kind of weird ordering for the numbering of roots, which we need to
// write back into depth-first ordering.
SetColumns(0);
SetQuitOnError(true);

procedure PrintReflectionTable(W)
    refltab := ReflectionTable(W);
    BFSOrder := {@ i : i in [1 .. Rank(W)] @};
    newtab := ZeroMatrix(Integers(), Rank(W), #refltab[1]);
    pos := 1;
    while pos le #BFSOrder do
        root := BFSOrder[pos];
        pos +:= 1;

        for s := 1 to Rank(W) do
            sroot := refltab[s][root];
            if sroot gt 0 then
                Include(~BFSOrder, sroot);
            end if;
            newtab[s, Index(BFSOrder, root)] := sroot gt 0 select Index(BFSOrder, sroot) else 0;
        end for;
    end while;
    CoxeterMatrix(W);
    printf "[%o]\n", newtab;
end procedure;

PrintReflectionTable(CoxeterGroup(GrpFPCox, "H3"));
PrintReflectionTable(CoxeterGroup(GrpFPCox, "H4"));

mat_353 := Matrix(Integers(), [
    [1, 3, 2, 2],
    [3, 1, 5, 2],
    [2, 5, 1, 3],
    [2, 2, 3, 1]
]);
PrintReflectionTable(CoxeterGroup(GrpFPCox, mat_353));

mat_rand := Matrix(Integers(), [
    [1, 3, 4, 2, 2, 2, 2],
    [0, 1, 3, 3, 2, 2, 2],
    [0, 0, 1, 9, 2, 2, 2],
    [0, 0, 0, 1, 8, 3, 2],
    [0, 0, 0, 0, 1, 3, 2],
    [0, 0, 0, 0, 0, 1, 5],
    [0, 0, 0, 0, 0, 0, 1]
]);
for i := 1 to Nrows(mat_rand) do for j := 1 to i-1 do mat_rand[i, j] := mat_rand[j, i]; end for; end for;
PrintReflectionTable(CoxeterGroup(GrpFPCox, mat_rand));



quit;

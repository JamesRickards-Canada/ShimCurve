intrinsic X0DNdata(DBound::RngIntElt, NBound::RngIntElt) -> Any
    {Outputs data table for the shimura curves X0(D;N), where D a quaternion discriminant and N coprime to D in the box D =< DBound and N <= NBound.}

    filename := Sprintf("./data/genera-tables/SignatureTableX0DN_%o_%o.txt", DBound, NBound);
    fprintf filename, "Glabel?all_degree1_points_known?autmuO_norms?bad_primes?cm_discriminants?coarse_class?coarse_class_num?coarse_index?coarse_label?coarse_num?conductor?curve_label?deg_mu?dims?discB?discO?fine_label?fine_num?fuchsian_index?galEnd?generators?genus?genus_minus_rank?gerbiness?has_obstruction?index?is_coarse?is_split?label?lattice_labels?lattice_x?level?level_is_squarefree?level_radical?log_conductor?models?mu_label?mults?name?newforms?nu2?nu3?nu4?nu6?num_bad_primes?num_known_degree1_noncm_points?num_known_degree1_points?obstructions?order_label?parents?parents_conj?pointless?power?psl2label?q_gonality?q_gonality_bounds?qbar_gonality?qbar_gonality_bounds?ram_data_elts?rank?reductions?scalar_label?simple?squarefree?torsion?trace_hash?traces \n";
    fprintf filename, "text?boolean?integer[]?integer[]?integer[]?text?integer?integer?text?integer?integer[]?text?integer?integer[]?integer?integer?text?integer?integer?text?integer[]?integer?integer?integer?smallint?integer?boolean?boolean?text?text[]?integer[]?integer?boolean?integer?numeric?smallint?text?integer[]?text?text[]?integer?integer?integer?integer?integer?integer?integer?integer[]?text?text[]?integer[]?boolean?boolean?text?integer?integer[]?integer?integer[]?numeric[]?integer?text[]?text?boolean?boolean?integer[]?bigint?integer[] \n\n";
    // we only want D square-free, with an even number of prime factors.
    for D in [D : D in [6..DBound] | MoebiusMu(D) eq 1] do
        // we want N that are coprime to D and square-free.
        for N in [N : N in [1..NBound] | GCD(D,N) eq 1] do
            Gamma := FuchsianGroup(QuaternionOrder(D,N) : VerifyEichler := false);
            // e := EllipticInvariants(Gamma);
            Glabel:="1.1";
            all_degree1_points_known:="\\N";
            autmuO_norms:="\\N";
            bad_primes:="\\N";
            cm_discriminants:="\\N";
            coarse_class:="a";
            coarse_class_num:="\\N";;
            coarse_index:="\\N";
            coarse_label:="\\N";
            coarse_num:=1;
            conductor:="\\N";
            curve_label:="\\N";
            deg_mu:=1;
            dims:="\\N";
            discB:=D;
            discO:=D*N;
            fine_label:="\\N";
            fine_num:="\\N";
            fuchsian_index:="\\N"; // the Index(Gamma) command doesn't work right away.. 
            galEnd:="\\N";
            generators:="\\N";
            genus:=Genus(Gamma);
            genus_minus_rank:="\\N";
            gerbiness:="\\N";
            has_obstruction:="\\N";
            index:="\\N";
            is_coarse:="\\N";
            is_split:="\\N";
            lattice_labels:="\\N";
            lattice_x:="\\N";
            level:=N;
            level_is_squarefree:= "T";
            level_radical:="\\N";
            log_conductor:="\\N";
            models:="\\N";
            mults:="\\N";
            name:="X_0(" cat Sprintf("%o,%o)", D,N);
            newforms:="\\N";
            nu2:=EllipticInvariants(Gamma)[1][2];
            nu3:=EllipticInvariants(Gamma)[2][2];
            nu4:=0;
            nu6:=0;
            num_bad_primes:="\\N";
            num_known_degree1_noncm_points:="\\N";
            num_known_degree1_points:="\\N";
            obstructions:="\\N";
            if N eq 1 then
                order_label:=Sprintf("%o",D);
            else
                order_label:=Sprintf("%o.%o",D,N);
            end if;
            mu_label:=order_label cat ".1";
            label:=mu_label cat coarse_label;
            parents:="\\N";
            parents_conj:="\\N";
            pointless:="\\N";
            power:="\\N";
            psl2label:="\\N";
            gonality_temp:=GonalityBoundListX0DN(D,N);
            q_gonality:=gonality_temp[1];
            q_gonality_bounds:=Sprintf("{%o,%o}",gonality_temp[2][1],gonality_temp[2][2]);
            qbar_gonality:=gonality_temp[3];
            qbar_gonality_bounds:=Sprintf("{%o,%o}",gonality_temp[4][1],gonality_temp[4][2]);
            ram_data_elts:="\\N";
            rank:="\\N";
            reductions:="\\N";
            scalar_label:="\\N";
            simple:="\\N";
            squarefree:="\\N";
            torsion:="\\N";
            trace_hash:="\\N";
            traces:="\\N";
            fprintf filename, Sprintf("%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o?%o \n",Glabel,all_degree1_points_known,autmuO_norms,bad_primes,cm_discriminants,coarse_class,coarse_class_num,coarse_index,coarse_label,coarse_num,conductor,curve_label,deg_mu,dims,discB,discO,fine_label,fine_num,fuchsian_index,galEnd,generators,genus,genus_minus_rank,gerbiness,has_obstruction,index,is_coarse,is_split,label,lattice_labels,lattice_x,level,level_is_squarefree,level_radical,log_conductor,models,mu_label,mults,name,newforms,nu2,nu3,nu4,nu6,num_bad_primes,num_known_degree1_noncm_points,num_known_degree1_points,obstructions,order_label,parents,parents_conj,pointless,power,psl2label,q_gonality,q_gonality_bounds,qbar_gonality,qbar_gonality_bounds,ram_data_elts,rank,reductions,scalar_label,simple,squarefree,torsion,trace_hash,traces);
        end for;
    end for;
    return Sprint("Table produced :)");
end intrinsic;
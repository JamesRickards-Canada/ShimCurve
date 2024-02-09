
intrinsic NormalizerPlusGenerators(O::AlgQuatOrd) -> SeqEnum 
  {return generators of the positive norm elements which normalize O}
  if Discriminant(O) eq 6 then 
    B6<i6,j6>:=QuaternionAlgebra<Rationals() | -1,3 >;
    B:=QuaternionAlgebra(O);
    tr,map:=IsIsomorphic(B6,B : Isomorphism:=true);
    assert tr;
    B6elliptic_elts:=[3*i6 + i6*j6, 1+i6, 3+3*i6+j6+i6*j6];
    Oelliptic_elts:=[ O!map(a) : a in B6elliptic_elts ];
    assert Set([ Norm(a) : a in Oelliptic_elts ]) eq {2,6,12};

    e2,e4,e6:=Explode(Oelliptic_elts);
    assert IsScalar(e6^6); assert IsScalar(e4^4); assert IsScalar(e2^2);
    assert IsScalar(&*Oelliptic_elts);
    return Oelliptic_elts;
  elif Discriminant(O) eq 10 then 
    //Elkies 
    B10<b,e>:=QuaternionAlgebra<Rationals() | -2,5 >;
    s2:=b;
    s2p:=2*e+5*b-b*e;
    s2pp:=5*b-b*e;
    s3:=2*b-e-1;

    B:=QuaternionAlgebra(O);
    tr,map:=IsIsomorphic(B10,B : Isomorphism:=true);
    assert tr;
    B10elliptic_elts:=[ s2,s2p,s2pp,s3];
    assert IsScalar(&*B10elliptic_elts);
    assert IsScalar(s2^2); assert IsScalar(s2p^2); assert IsScalar(s2pp^2); assert IsScalar(s3^3);
    Oelliptic_elts:=[ O!map(a) : a in B10elliptic_elts ];
    //assert Set([ Norm(a) : a in Oelliptic_elts ]) eq {2,6,12};
    return Oelliptic_elts;
  elif Discriminant(O) eq 15 then 
    B15<c,e>:=QuaternionAlgebra<Rationals() | -3,5 >;
    s2:=4*c-3*e;
    s2p:=5*c-3*e-c*e;
    s2pp:=20*c-9*e-7*c*e;
    s6:=3+c;

    B:=QuaternionAlgebra(O);
    tr,map:=IsIsomorphic(B15,B : Isomorphism:=true);
    assert tr;
    B15elliptic_elts:=[ s2,s2p,s2pp,s6 ];
    assert IsScalar(&*B15elliptic_elts);
    assert IsScalar(s2^2); assert IsScalar(s2p^2); assert IsScalar(s2pp^2); assert IsScalar(s6^6);

    Oelliptic_elts:=[ O!map(a) : a in B15elliptic_elts ];
    //assert Set([ Norm(a) : a in Oelliptic_elts ]) eq {2,6,12};
    return Oelliptic_elts;

  else
    return "oops, not written for this discriminant yet";
  end if;
end intrinsic;




intrinsic SemidirectToNormalizer(O::AlgQuatOrd,mu::AlgQuatOrdElt,h::AlgQuatEnhElt) -> AlgQuatProjElt
  {the map from the semidirect product to the normalizer.}
  w:=(h`element[1])`element;
  x:=h`element[2];
  return Parent(h`element[1])!(w*x);
end intrinsic;



intrinsic SemidirectToNormalizerKernel(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> SeqEnum 
  {return the kernel of the map form the enhanced semidirect product to N_B^x(O). 
  It is necessarily cyclic and the second value is the generator of the group}
  B:=QuaternionAlgebra(O);
  Ocirc:=EnhancedSemidirectProduct(O);
  AutFull, autmuOseq := Aut(O,mu);
  Oxcyc:= [ (1/Integers()!Sqrt(Norm(a`element)))*a`element : a in autmuOseq | IsSquare(Norm(a`element)) ];
  ker:=[ Ocirc!<x,x^-1> : x in Oxcyc ];
  assert #ker in [1,2,3];
  assert Set([ Norm(e) eq 1 : e in Oxcyc ]) eq Set([true]);
  if #ker eq 1 then 
    assert ker[1] eq Ocirc!<B!1,O!1> or ker[1] eq Ocirc!<B!1,-O!1>;
    return [ Ocirc!<B!1,O!1>,Ocirc!<B!1,-O!1> ],Ocirc!<B!1,-O!1>;
  else 
    gen:=[ e : e in ker | Order(e) eq 2*#ker ];
    assert #gen eq 1;
    gen:=gen[1];
    newker:=[ gen^i : i in [1..Order(gen)] ];
    assert #Set(newker) eq Order(gen);
    //assert its cyclic in GL4
    return newker,gen;
  end if;
end intrinsic;

intrinsic SemidirectToNormalizerKernel(O::AlgQuatOrd,mu::AlgQuatElt) -> SeqEnum 
  {return the kernel of the map form the enhanced semidirect product to N_B^x(O). 
  It is necessarily cyclic and the second value is the generator of the group}
  return SemidirectToNormalizerKernel(O,O!mu);
end intrinsic;

intrinsic NormalizerToAutmuO(O::AlgQuatOrd,mu::AlgQuatOrdElt,a::AlgQuatOrdElt) -> AlgQuatEnhElt 
  {Lift an element a of the Normalizer of O to the enhanced semidirect product, which is well defined up to 
  the kernel of this map (given by SemidirectToNormalizerKernel)}
  Ocirc:=EnhancedSemidirectProduct(O);
  AutFull,autmuOseq:=Aut(O,mu);
  ker,kergen:=SemidirectToNormalizerKernel(O,mu);


  B:=QuaternionAlgebra(O);
  BxmodQx:=QuaternionAlgebraModuloScalars(B);
  proja:=BxmodQx!(B!a);
  orda:=Order(proja);

  //[ elt : elt in autmuOseq | elt in ker ];

  assert a^2/Norm(a) in O;
  assert Norm(a) gt 0;
  W:=[];
  for w in autmuOseq do 
    if IsSquare(Rationals()!Abs(Norm((w`element)^-1*a))) then
      tr,c:=IsSquare(Rationals()!Abs(Norm((w`element)^-1*a)));
      x:=(1/c)*((w`element)^-1)*a;
      assert x in O;
      assert Norm(x) in {1,-1};
      ell:=Ocirc!<w,O!x>;
      if Min([ i : i in [1..orda] | ell^i in ker]) eq orda then 
        Append(~W,ell);
        //return ell;
      end if;
    end if;
  end for;
  return W[1];
end intrinsic;

intrinsic NormalizerToAutmuO(O::AlgQuatOrd,mu::AlgQuatElt,a::AlgQuatElt) -> AlgQuatEnhElt 
  {Lift an element a of the Normalizer of O to the enhanced semidirect product, which is well defined up to 
  the kernel of this map (given by SemidirectToNormalizerKernel)}
  return NormalizerToAutmuO(O,O!mu,O!a);
end intrinsic;


intrinsic NormalizerToAutmuO(O::AlgQuatOrd,mu::AlgQuatElt,a::AlgQuatOrdElt) -> AlgQuatEnhElt 
  {Lift an element a of the Normalizer of O to the enhanced semidirect product, which is well defined up to 
  the kernel of this map (given by SemidirectToNormalizerKernel)}
  return NormalizerToAutmuO(O,O!mu,a);
end intrinsic;

intrinsic NormalizerToAutmuO(O::AlgQuatOrd,mu::AlgQuatOrdElt,a::AlgQuatElt) -> AlgQuatEnhElt 
  {Lift an element a of the Normalizer of O to the enhanced semidirect product, which is well defined up to 
  the kernel of this map (given by SemidirectToNormalizerKernel)}
  return NormalizerToAutmuO(O,mu,O!a);
end intrinsic;




intrinsic NormalizerPlusGeneratorsEnhanced(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> Tup 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  ker,kergen:=SemidirectToNormalizerKernel(O,mu);
  Ocirc:=EnhancedSemidirectProduct(O);
  Nplus:=NormalizerPlusGenerators(O);
  return [ Ocirc!NormalizerToAutmuO(O,O!mu,O!a) : a in NormalizerPlusGenerators(O) ] cat [Ocirc!kergen];
end intrinsic;

intrinsic NormalizerPlusGeneratorsEnhanced(O::AlgQuatOrd,mu::AlgQuatElt) -> Tup 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  return NormalizerPlusGeneratorsEnhanced(O,O!mu);
end intrinsic;

intrinsic NormalizerPlusGeneratorsEnhanced(O::AlgQuatOrd,del::RngIntElt) -> Tup 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  tr,mu:=HasPolarizedElementOfDegree(O,del);
  return NormalizerPlusGeneratorsEnhanced(O,O!mu);
end intrinsic;



intrinsic NormalizerPlusGeneratorsGL4modN(O::AlgQuatOrd,mu::AlgQuatOrdElt,N::RngIntElt) -> SeqEnum 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  return [ EnhancedElementInGL4modN(g,N) : g in NormalizerPlusGeneratorsEnhanced(O,mu) ];
end intrinsic;

intrinsic NormalizerPlusGeneratorsGL4modN(O::AlgQuatOrd,mu::AlgQuatElt,N::RngIntElt) -> SeqEnum 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  return [ EnhancedElementInGL4modN(g,N) : g in NormalizerPlusGeneratorsEnhanced(O,O!mu) ];
end intrinsic;


intrinsic NormalizerPlusGeneratorsGL4modN(O::AlgQuatOrd,del::RngIntElt,N::RngIntElt) -> SeqEnum 
  {return generators of the positive norm elements which normalize O in the enhanced semidirect product}
  tr,mu:=HasPolarizedElementOfDegree(O,del);
  return [ EnhancedElementInGL4modN(g,N) : g in NormalizerPlusGeneratorsEnhanced(O,mu) ];
end intrinsic;





intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatOrdElt) -> SeqEnum 
  {return the elliptic elements}
  Ocirc:=EnhancedSemidirectProduct(O);
  return [ Ocirc!NormalizerToAutmuO(O,mu,a) : a in NormalizerPlusGenerators(O) ];
end intrinsic;

intrinsic EnhancedEllipticElements(O::AlgQuatOrd,mu::AlgQuatElt) -> SeqEnum
  {return the elliptic elements of the enhanced semidirect product}

  return EnhancedEllipticElements(O,O!mu);
end intrinsic;






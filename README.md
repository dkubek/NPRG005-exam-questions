# NeprocProg Cure

[TOC]

## Prolog

### Topologické uspořádání grafu

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095&sid=fe143536d7d0b2e925781a412fafbdc8)

Je dán orientovaný graf G pomocí seznamů sousedů. Zjistěte, jestli lze graf G topologicky uspořádat a pokud ano, vydejte seznam vrcholů v topologickém pořadí.

Příklad:
```prolog
?- topo([a-[],b-[a,c],c-[a],d-[a,c]],Usp).
Usp=[b,d,c,a]
```

 1) Definujte příslušný predikát topo/2 v jazyce Prolog.
 2) Odhadněte časovou složitost vašeho řešení. Odhad zdůvodněte.
 3) Jsou některé z vašich predikátů koncově rekurzivní ? Pokud ano, vysvětlete, které to jsou, a jaký to má význam. Pokud ne, vysvětlete, zdali by se dal některý takto upravit. 

Řešení:
[TODO]

### Diskrepanční vrstvy

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095&sid=fe143536d7d0b2e925781a412fafbdc8)

Napište predikát ``diskr/2``, který dostane binární strom (s konstruktory ``t/3`` a ``nil/0``) a vrátí seznam seznamů vrcholů stromu, kde v jednom vnitřním seznamu jsou všechny vrcholy, ke kterým se při průchodu od kořene dostaneme se stejným počtem kroků doprava. Vnější seznam je od nejlevější vrstvy, na pořadí ve vnitřních seznamech nezáleží.

Příklad:
```
?- diskr(t( t(t(nil,a,nil),b,t(nil,c,nil)),
            d,
            t(t(nil,e,t(nil,f,nil)),
              g,
              t(nil,h,t(nil,i,nil)) )), V).
V = [[a,b,d],[c,g,e],[f,h],[i]]
```

1. Definujte příslušný predikát ``diskr/2.``
2. Je ve vašem řešení použit řez (!) nebo negace? Pokud ano, změní se něco, když řez / negaci vypustíme? Pokud ne, dal by se řez / negace někde smysluplně využít?
3. Lze u predikátu ``diskr/2`` obrátit směr výpočtu? Podrobněji: dle příkladu předpokládáme volání diskr(+,-). Bude fungovat i volání diskr(-, +), tj. zadáme seznam diskrepančních vrstev, a na výstupu obdržíme strom? Vysvětlete. 

Řešení:
```prolog=
diskr(Tree, V) :-
    diskr_(Tree, 0, NodeRightCount),
    collect(NodeRightCount, 0, V).

diskr_(nil, _, []).
diskr_(t(Left, Node, Right), RightCount, NodeRightCount) :-
    diskr_(Left, RightCount, LAns),
    NewCount is RightCount + 1,
    diskr_(Right, NewCount, RAns),
    append(LAns, [Node-RightCount | RAns], NodeRightCount).


is_count(Count, _-Count).
get_value(Value-_, Value).

collect([], _, []).
collect(Pairs, Count, Ans) :-
    include(is_count(Count), Pairs, RightCount),
    (
        RightCount \= []
    ->
        maplist(get_value, RightCount, Vals),
        NewCount is Count + 1,
        collect(Pairs, NewCount, TmpAns),
        Ans = [Vals | TmpAns]
    ;
        Ans = []
    ).


% test data
test_tree(
    t(
        t(
            t(nil, a, nil),
            b,
            t(nil, c, nil)
        ),
        d,
        t(
            t(
                nil,
                e,
                t(nil, f, nil)
            ),
            g,
            t(
                nil,
                h,
                t(nil, i, nil)
            )
        )
    )
).

```

### Generování binárních stromů

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089&sid=fe143536d7d0b2e925781a412fafbdc8)

Cílem úlohy je definovat predikát allTrees, který pro daný seznam hladin vygeneruje všechny možné binární stromy.

 - Hladinou rozumíme seznam prvků, které se nacházejí ve stejné hloubce
 - Můžete předpokládat, že každá hladina má nanejvýš dvojnásobek prvků předchozí hladiny (ale může jich mít méně).
 - Hladiny vygenerovaného stromu musejí odpovídat hladinám specifikovaných ve vstupním seznamu.

Např. pro seznam ``[[1],[2,3],[4]]`` dostaneme následující 4 stromy:

```
   1
 2   3
4

   1
 2   3
  4

   1
 2   3
    4

   1
 2   3
      4
```

 1. Popište zvolenou reprezentaci binárních stromů.
 2. Definujte predikát ``allTrees/2``.
 3. Stručně vysvětlete, proč je vaše definice korektní.
 4. Lze vaší definici použít opačným směrem? Tj. nalezne váš predikát seznam hladin pokud specifikujete pouze výsledný strom? Vysvětlete.

Řešení:
```prolog=
level_to_forest([], [], []).
level_to_forest([X | Xs], [Left, Right | Rest], [t( Left, X, Right ) | Ans]) :-
    level_to_forest( Xs, Rest, Ans ).
level_to_forest([X | Xs], [Left | Rest],        [t( Left, X, nil ) | Ans]) :-
    level_to_forest( Xs, Rest, Ans ).
level_to_forest([X | Xs], [Right | Rest],       [t( nil, X, Right ) | Ans]) :-
    level_to_forest( Xs, Rest, Ans ).
level_to_forest([X | Xs], Trees,                [t( nil, X, nil ) | Ans]) :-
    level_to_forest( Xs, Trees, Ans ).

all_trees(Levels, SingleTree) :-
    reverse(Levels, ReversedLevels),
    Forest = [],

    all_trees(ReversedLevels, Forest, [SingleTree]).

all_trees([], Forest, Forest).
all_trees([Level | Levels], Forest, Ans) :-
    level_to_forest(Level, Forest, NewForest),
    all_trees(Levels, NewForest, Ans).
```


### Bipartitní rozklad grafu

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089&sid=fe143536d7d0b2e925781a412fafbdc8)

Je zadán neorientovaný graf *G* a množina vrcholů *M*. Zjistěte, zda *M* a doplněk *M* tvoří bipartitní rozklad grafu *G* (tj. každá hrana grafu má právě jeden koncový vrchol v množině *M*). Pokud ano, vydejte druhou množinu rozkladu.

```prolog
?- bip([a-[c,d], b-[d], c-[a], d-[a,b]], [a,b], D).
    D = [c,d]

?- bip([a-[c,d], b-[d], c-[a], d-[a,b]], [b,c], D).
    false
```

 1. Definujte predikát ``bip/3``.
 2. Napište o jednotlivých predikátech ve vašem řešení, zda jsou koncově rekurzivní.

Řešení:
```prolog=
% collect_nodes(+Graph, -Nodes) is true when Nodes are all the nodes in Graph
% in sorted order.
collect_nodes(Graph, Nodes) :-
    collect_nodes(Graph, [], NodesDup),
    sort(NodesDup, Nodes),
    !.

collect_nodes([], Acc, Acc).
collect_nodes([Node-Neighbours | Ns], Acc, Ans) :-
    append([Node | Neighbours], Acc, NewAcc),
    collect_nodes(Ns, NewAcc, Ans).


% difference(+List1, +List2, -List3) is true if List3 contains all nodes of List1
% except the elements in List2
difference(List, [], List).
difference(List, [X | Xs], Out) :-
    exclude(=(X), List, Tmp),
    difference(Tmp, Xs, Out),
    !.


bip(Graph, V, U) :-
    collect_nodes(Graph, Nodes),
    difference(Nodes, V, U),
    is_bipartite(Graph, V, U).

is_bipartite([], _, _).
is_bipartite([Node-Neighbours | NNs], V, U) :-
    (
        member(Node, V)
    ->
        Partition = V
    ;
        Partition = U
    ),
    member(Node, Partition),
    maplist(does_not_contain(Partition), Neighbours),
    is_bipartite(NNs, V, U),
    !.

does_not_contain(List, X) :- \+ member(X, List).
```

### Problém truhláře

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078&sid=fe143536d7d0b2e925781a412fafbdc8)

Truhlář má dostatek trámů délky ``D`` a seznam ``Xs`` délek trámů, které potřebuje nařezat. V seznamu ``Xs`` se délky mohou opakovat.

Cílem problému je sestavit predikát ``rezy(+D, +Xs, -N, -Vss)``, který

 - rozdělí požadované délky do skupin, které se mají nařezat z jednoho trámu
 - truhlář přitom používá hladový algoritmus, tj. pro každou délku použije první trám, z něhož lze ještě požadovanou délku odřezat
 - vrátí celkový počet řezaných trámů N
 - a seznam seznamů Vss (délky N), jehož každý prvek reprezentuje dělení jednoho trámu (případný zbytek se neuvádí).

```prolog
?- rezy(5,[3,2,2,2,2,1,4], N, V).
N=4, V=[[3,2],[2,2,1],[2],[4]]
```

1. Definujte predikát ``rezy/4.`` Definice případných pomocných predikátů prosím opatřete vysvětlujícím komentářem.
2. Je některý z vašich predikátů koncově rekurzivní? Pokud ano, vysvětlete, který to je a jaký to má význam.
3. Pokud ne, dal by se některý takto upravit? Odpověď prosím zdůvodněte.
    
Řešení:
```prolog=
second(_-X, X).

rezy(Length, Xs, N, Vss) :-
    rezy_(Length, Xs, [], TRVss),
    maplist(second, TRVss, RVss),
    maplist(reverse, RVss, Vss),
    length(Vss, N),
    !.

rezy_(_, [], Acc, Acc).
rezy_(Length, [ X | Xs ], Acc, Ans) :-
    greedy_extend_cut(Length, Acc, X, NewAcc),
    rezy_(Length, Xs, NewAcc, Ans).

greedy_extend_cut(_, [], NextCut, [ NextCut-[NextCut] ]).
greedy_extend_cut(MaxLength, [ Total-Cuts | Cs ], NextCut, Ans) :-
    NextCut =< MaxLength,
    Free is MaxLength - Total,
    (
        NextCut =< Free
    ->
        NewTotal is Total + NextCut,
        Ans = [ NewTotal-[NextCut | Cuts ] | Cs ]
    ;
        greedy_extend_cut(MaxLength, Cs, NextCut, Ans_),
        Ans = [ Total-Cuts | Ans_]
    ),
    !.
```
    
### Systém různých reprezentantů

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078&sid=fe143536d7d0b2e925781a412fafbdc8)

Je zadán seznam množin ``Mss``. Chceme všemi možnými způsoby vybrat a vrátit v seznamu reprezentanty daných množin v odpovídajícím pořadí s podmínkou, že konkrétní reprezentanti v jednom výběru jsou různí.

Příklad:
```prolog
?- reprezentanti([[1],[1,2,3],[1,3,4]], R).
R = [[1,2,3],[1,2,4],[1,3,4]]
```

1. Sestavte predikát ``reprezentanti(+Mss, -Rss)``.
2. Stručně vysvětlete, proč je vaše definice korektní.
3. Je ve vašem programu použit řez ``(!)`` ? Jde o řez červený (mění deklarativní význam programu) či zelený (nemění d.v.)? Pokud ne, je řez nezbytný pro definici některého vestavěného predikátu / operátoru, který jste ve vašem řešení použili? Jde o řez červený (mění deklarativní význam programu) či zelený (nemění d.v.)?

```prolog=
cons(X, Xs, [X | Xs]).

extend_one(X, Ys, Ans) :-
    exclude(member(X), Ys, Tmp),
    maplist(cons(X), Tmp, Ans).

extend_with([], _, Acc, Acc).
extend_with([X | Xs], Ys, Acc, Ans) :-
    extend_one(X, Ys, Tmp),
    append(Tmp, Acc, NewAcc),
    extend_with(Xs, Ys, NewAcc, Ans).


reprezentanti(Xss, Ans) :-
    reprezentanti(Xss, [[]], Ans).

reprezentanti([], Acc, Acc).
reprezentanti([Xs | Xss], Acc, Ans) :-
    extend_with(Xs, Acc, [], NewAcc),
    reprezentanti(Xss, NewAcc, Ans).
```

###  Hammerstein

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066&sid=fe143536d7d0b2e925781a412fafbdc7)

Profesor Hammerstein definoval predikat ``setrid/2`` takto:

```prolog
% setrid(+Xs,-Ys) :- Ys je seznam přirozených čísel ze seznamu Xs setříděný
% vzestupně
setrid(Xs,Ys) :-
    append(A,[H1,H2|B],Xs),
    H1 > H2,
    !,
    append(A,[H2,H1|B],Xs1),
    setrid(Xs1,Ys).
```

zapomněl však na klauzuli, která definuje bázi rekurze.

1. Doplňte jednu (opravdu jen jednu) chybějící klauzuli za uvedené pravidlo tak, aby výsledná procedura korektně setřídila vstupní seznam přirozených čísel. Na výstupu bychom měli obdržet jen jediné řešení.
2.  V definici pravidla je použit řez (!). Jde o zelený (nemění deklarativní význam) či červený řez (mění d.v.) ? Vysvětlete! Obsahuje některá z vašich klauzulí, (doplněná v(a) nebo (b)) zelený či červený řez?
3. Jaký známý třídící algoritmus výše uvedený kód implementuje? Pokud neznáte název, můžete alespoň slovně popsat, jak setrid/2 funguje.
4. *VOLITELNE*: Lze u procedury ``setrid/2`` obrátit směr výpočtu?
```prolog
setrid(-Xs,+Ys) :- Xs je seznam přirozených čísel ze seznamu Ys setříděný vzestupně
```

Pokud ne, šel by kód jednoduše upravit tak, aby se výsledný predikát (pojmenovaný třeba ``setrid2/2``) dal korektně volat oběma způsoby?

Řešení:
```prolog=
% setrid(+Xs,-Ys) :- Ys je seznam přirozených čísel ze seznamu Xs setříděný
% vzestupně
% Bubble sort
setrid(Xs, Ys) :-
    append(A, [H1,H2|B], Xs),
    H1 > H2,
    !, % cerveny
    append(A, [H2,H1|B], Xs1),
    setrid(Xs1, Ys).
setrid(Xs, Xs).
```

### Cestovatel

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066&sid=fe143536d7d0b2e925781a412fafbdc7)

Do země Mobilia, v níž je každý občan vybaven chytrým telefonem, přicestoval Cestovatel, nakažený virovým onemocněním. Všichni ostatní byli přitom ještě zdraví. Můžeme předpokládat, že virus se přenese z jedné osoby na druhou, pokud spolu strávili ve vzdálenosti menší než 2m alespoň čas ``K``, kde ``K`` je známá kritická hodnota. Díky chytrým telefonům máme pro každého občana Mobilie seznam záznamů jeho kontaktů, kde každý takový záznam pro osobu ``A`` obsahuje identifikaci osoby ``B``, která se k němu přiblížila do vzdálenosti ``< 2m`` čas setkání a délku setkání.

Cílem je sestavit program, který na základě takových záznamů vrátí seznam infikovaných osob.

1. V jazyce Prolog popište datovou strukturu pro reprezentaci jednoho záznamu kontaktu občana Mobilie popsaného výše.
2. V jazyce Prolog navrhněte reprezentaci položek VstupníhoSeznamu, přičemž každá položka bude obsahovat indentifikaci občana Mobilie a seznam záznamů jeho kontaktů.
3. Sestavte predikát ``inf/4``, který obdrží
```
VstupníSeznam
identifikaci Cestovatele
kritickou hodnotu K
```
a vrátí seznam infikovaných.

U každého pomocného predikátu prosím v poznámce popište jeho význam.

*Volitelné:* výstupní seznam můžete uspořádat dle délky kontaktu s infikovanými do nerostoucí posloupnosti.

4. Odhadněte časovou složitost vašeho řešení.
5. Je některý z vašich predikátů koncově rekurzivní ? Pokud ano, vysvětlete, jaký to má význam. Pokud ne , dal by se některý takto upravit?

Řešení:
```prolog=
% contact(ID, Time, Length)

% [ID1-[ contact(...), ... ], ID2-[ contact(...), ... ], ...]

lookup_contact([ID-Contacts | _], ID, Contacts) :- !.
lookup_contact([ _ | Cs ], ID, Contacts) :-
    lookup_contact(Cs, ID, Contacts).

was_infected(TimeInfected, K, contact(_, Time, Length)) :-
    Time >= TimeInfected,
    Length >= K.

get_id_time_pair(contact(ID, Time), ID-Time).

get_infected(ContactsList, K, ID, TimeInfected, Infected) :-
    lookup_contact(ContactsList, ID, Contacts),
    include(was_infected(TimeInfected, K), Contacts, InfectedCs),
    maplist(get_id_time_pair, InfectedCs, Infected).

contained_in(InfectedIDs, ID-_) :- member(ID, InfectedIDs).

first(X-_, X).

inf(ContactsList, TravellerID, K, Infected) :-
    TimeInfected = 0,
    ToProcess = [ TravellerID-TimeInfected ],
    Acc = [ ],
    inf(ContactsList, K, ToProcess, Acc, Infected).

inf(_, _, [], Acc, Acc).
inf(ContactsList, K, [ ID-TimeInfected | Ps ], Acc, Ans) :-
    get_infected(ContactsList, K, ID, TimeInfected, Infected),

    exclude(contained_in(Acc), Infected, NewInfected),
    append(Ps, NewInfected, NewToProcess),

    maplist(first, NewInfected, NewIDs),
    append(NewIDs, Acc, NewAcc),

    inf(ContactsList, K, NewToProcess, NewAcc, Ans).
```


### Generování hodnot výrokových proměnných

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977&sid=fe143536d7d0b2e925781a412fafbdc8)

Definujte binární predikát ``aspon2/2``, který
- obdrží seznam výrokových proměnných (reprezentovaných atomy), v němž je každá proměnná ohodnocena hodnotou true nebo false
- vrátí seznam všech takových ohodnocení týchž proměnných, v němž se každé ohodnocení bude od vstupního lišit v hodnotách alespoň 2 proměnných.

Příklad:
```prolog
?- aspon2([x1-true, x2-false, y-true], V).
  V =  [  [x1-false, x2-true, y-true],
          [x1-false, x2-false, y-false],
          [x1-true, x2-true, y-false],
          [x1-false, x2-true, y-false] ]
```

Řešení:
```prolog=
cons(X, Xs, [X | Xs]).

% diff(Modified, Original, N).
diff_ord([], [], 0).
diff_ord([_-V | T1], [_-not(V) | T2], AnsN) :-
    diff_ord(T1, T2, N),
    AnsN is N + 1.
diff_ord([_-V | T1], [_-V | T2], N) :-
    diff_ord(T1, T2, N).

diff_less_then_2(Original, Modified) :-
    diff_ord(Original, Modified, N),
    N < 2.

aspon_2(Values, X) :-
    subset_change(Values, All),
    exclude(diff_less_then_2(Values), All, X).

subset_change([], [[]]).
subset_change([Var-Value | Values], Ans) :-
    subset_change(Values, Tmp),

    maplist(cons(Var-Value), Tmp, NotChanged),
    maplist(cons(Var-not(Value)), Tmp, Changed),

    append(NotChanged, Changed, Ans).
```

### Trojúhelníky v grafu

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977&sid=fe143536d7d0b2e925781a412fafbdc8)

Graf je zadán jako seznam svých vrcholů se seznamy sousedů (viz příklad). Definujte binární predikát ``troj(+Graf, -SeznamTrojuhelniku)`` který k takovému grafu vrátí seznam všech jeho trojúhelníků. Ve výsledném seznamu by se každý trojúhelník měl vyskytovat právě jednou (``t(a,b,c)``, ``t(b,c,a)`` a ``t(c,a,b)`` jsou stejné trojúhelníky).

Příklad:
```prolog
?- troj([a-[b,c,d],b-[a,c],c-[a,b,d],d-[a,c],e-[]], S).
     S = [t(a,b,c), t(a,c,d)]
```
     
Řešení:
```prolog=
is_edge(Graph, From, To) :-
    member(From-Neighbours, Graph),
    member(To, Neighbours).

neighbours(Graph, Node, Neighbours) :-
    member(Node-Neighbours, Graph).

first(X-_, X).
lift_list(X, [X]).

all_paths(_, 0, []).
all_paths(Graph, N, Paths) :-
    maplist(first, Graph, Nodes),
    maplist(lift_list, Nodes, Acc),
    all_paths(Graph, N, 1, Acc, Paths).

all_paths(_, N, N, Acc, Acc) :- !.
all_paths(Graph, MaxLen, Len, Paths, Ans) :-
    extend_paths(Graph, Paths, NewAcc),
    NewLen is Len + 1,
    all_paths(Graph, MaxLen, NewLen, NewAcc, Ans),
    !.


extend_paths(Graph, Paths, Ans) :-
    extend_paths(Graph, Paths, [], Ans).

push(Xs, X, [X | Xs]).

extend_paths(_, [], Acc, Acc).
extend_paths(Graph, [ Path | Ps ], Acc, Ans) :-
    Path = [ N | _ ],
    neighbours(Graph, N, Neighbours),
    maplist(push(Path), Neighbours, Tmp),

    append(Tmp, Acc, NewAcc),

    extend_paths(Graph, Ps, NewAcc, Ans).

is_triangle([A, B, C, A]) :-
    A \= B,
    B \= C,
    A \= C.

to_triangle([A, B, C, A], t(A, B, C)).

is_congruent(t(A, B, C), t(A, B, C)).
is_congruent(t(A, B, C), t(A, C, B)).
is_congruent(t(A, B, C), t(B, A, C)).
is_congruent(t(A, B, C), t(B, C, A)).
is_congruent(t(A, B, C), t(C, A, B)).
is_congruent(t(A, B, C), t(C, B, A)).

deduplicate(Triangles, Ans) :-
    deduplicate(Triangles, [], Ans).

deduplicate([], Acc, Acc).
deduplicate([ Triangle | Ts ], Acc, Ans) :-
    exclude(is_congruent(Triangle), Ts, Filtered),
    deduplicate(Filtered, [Triangle | Acc], Ans).

troj(Graph, Triangles) :-
    all_paths(Graph, 4, Paths),
    include(is_triangle, Paths, Filtered),
    maplist(to_triangle, Filtered, Tmp),
    deduplicate(Tmp, Triangles).
```


### Generování výrokových formulí

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969&sid=fe143536d7d0b2e925781a412fafbdc8)

Formule výrokového počtu jsou sestavené z (výrokových) proměnných ve funktoru ``var/1`` a logických spojek negace, konjunkce a disjunkce (bez konstant). Dále máte dány v argumentech predikátu ``gen/3`` číslo ``k`` pro velikost formule a seznam jmen proměnných. Generujte backtrackingem všechny logické formule (každou jednou), které obsahují proměnné ze seznamu a ve kterých je počet spojek a výskytů proměnných dohromady právě ``k``.

Definujte predikát ``gen(+K, +Jmena, -Fle)``. Na pořadí generovaných formulí nezáleží, ale měli byste vygenerovat každou právě jednou. K řešení není potřeba predikát ``=../2`` (univ).

Příklad:
```prolog
?- gen(4,[p],F).

F = not(not(not(var(p))));
F = not(and(var(p),var(p)));
F = not(or(var(p),var(p)));
F = and(not(var(p)),var(p));
F = and(var(p),not(var(p)));
F = or(not(var(p)),var(p));
F = or(var(p),not(var(p)));
false.
```

Řešení:
```prolog=
gen(K, Vars, F) :-
    length(Slots, K),
    gen_(Slots, Vars, F).

gen_([_], Vars, var(V)) :-
    member(V, Vars).
gen_([_ | Ss], Vars, not(F)) :-
    gen_(Ss, Vars, F).
gen_([_ | Ss], Vars, Ans) :-
    append(Left, Right, Ss),
    gen_(Left, Vars, F1),
    gen_(Right, Vars, F2),
    (
        Ans = and(F1, F2)
    ;
        Ans = or(F1, F2)
    ).
```

### Koncepty

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969&sid=fe143536d7d0b2e925781a412fafbdc8)

Jeden objekt je zadán uspořádaným seznamem dvojic klíč-hodnota. Na vstupu máte seznam objektů. Napište proceduru ``koncept/2``, která vyrobí nejmenší koncept zahrnující všechny vstupní objekty. Koncept je seznam dvojic ``klíč-seznam_hodnot``. Koncept zahrnuje objekt, pokud koncept má všechny klíče objektu a v seznamu hodnot příslušného klíče u konceptu je obsažena hodnota klíče u objektu. Pokud objekt nějaký klíč konceptu nemá, bude v seznamu hodnot konceptu hodnota ``nedef``.

Příklad:
```prolog
?- koncept([ [barva-modra, motor-diesel, pocet_kol-6],
             [barva-bila, motor-plyn, pocet_mist-40],
             [motor-elektro, pocet_mist-5] ],
             Koncept).
Koncept = [ barva-[modra,bila,nedef],
            motor-[diesel,plyn,elektro],
            pocet_kol-[6,nedef],
            pocet_mist-[40,5,nedef] ]
```

Řešení:
```prolog=
collect_attributes(Objects, Attributes) :-
    collect_attributes(Objects, [], Attributes).

collect_attributes([], Acc, Ans) :-
    sort(Acc, Ans).
collect_attributes([ Object | Os ], Acc, Ans) :-
    collect_attributes_one(Object, Attrs),
    append(Attrs, Acc, NewAcc),
    collect_attributes(Os, NewAcc, Ans).

collect_attributes_one(Object, Attributes) :-
    collect_attributes_one(Object, [], Attributes).

collect_attributes_one([], Acc, Acc).
collect_attributes_one([Key-_ | Ps], Acc, Ans) :-
    collect_attributes_one(Ps, [ Key | Acc ], Ans).

koncept(Objects, Concepts) :-
    collect_attributes(Objects, Attrs),
    koncept(Objects, Attrs, [], Concepts).

koncept([], _, Concepts, Concepts).
koncept([ Object | Os ], Attrs, Concepts, Ans ) :-
    extend_concepts(Attrs, Object, Concepts, NewConcepts),
    koncept( Os, Attrs, NewConcepts, Ans ).

extend_concepts([], _, Concepts, Concepts).
extend_concepts([ Attr | Attrs ], Objects, Concepts, Ans) :-
    (
        select(Attr-Value_, Objects, RestObjects)
    ->
        Value = Value_,
        NewObjects = RestObjects
    ;
        Value = nedef,
        NewObjects = Objects
    ),
    extend_concept(Attr-Value, Concepts, NewConcepts),
    extend_concepts(Attrs, NewObjects, NewConcepts, Ans).


extend_concept(Attr-Value, Concepts, NewConcepts) :-
    (
        select(Attr-Values, Concepts, RestConcepts)
    ->
        set_add(Values, Value, NewValues),
        NewConcept = Attr-NewValues,
        NewConcepts = [ NewConcept | RestConcepts ]
    ;
        NewValues = [ Value ],
        NewConcept = Attr-NewValues,
        NewConcepts = [ NewConcept | Concepts ]
    ).


set_add(List, Element, Out) :-
    (
        member(Element, List)
    ->
        Out = List
    ;
        Out = [ Element | List ]
    ).
```

### Překrytí segmentů

Zdroj: [MFF Forum: Zkouška 10.6.2019 (Dvořák + Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11954)

Máte dány dva řetězce, u kterých nevíte jejich vzájemnou orientaci. Najděte a vydejte v seznamu všechna jejich vzájemná neprázdná překrytí.
Příklad:

Kód: Vybrat vše

?- prekryti([a,t,c,t,c],[c,t,c,c], V).
V = [a,t,c,t,c,t,c,c],[a,t,c,t,c,c],[a,t,c,t,c,c,t,c]]

### Neporovnatelné prvky částečně uspořádané množiny

Zdroj: [MFF Forum: Zkouška 10.6.2019 (Dvořák + Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11954)

Částečně uspořádaná množina je popsána seznamem termů tvaru x -> y s významem x pokrývá y (tj. x > y a současně x ≥ z ≥ y implikuje x = z nebo y = z).
Definujte predikát nepor/2, který k takto zadané množině vrátí seznam všech dvojic vzájemně neporovnatelných prvků (tj. dvojic x,y takových, že neplatí x ≥ y ani x ≤ y).
Příklad:

Kód: Vybrat vše

?- nepor([a->b, a->c, b->d, e->f], N).
N = [a-e,a-f,b-c,b-e,b-f,c-d,c-e,c-f,d-e,d-f]

### Lexikograficky předchozí permutace

Zdroj: [MFF Forum: Zkouška 21.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11747)


Nalezněte lexikograficky předchozí permutaci:(Pokud není tak false)

Kód: Vybrat vše

prev([1,2,6,3,4,5,7],V).
   V = [1,2,5,7,6,4,3] 
   


## Haskell

### Největší součet souvislé podposloupnosti

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095&sid=fe143536d7d0b2e925781a412fafbdc8)

1. Pro zadanou posloupnost čísel najděte spojitý úsek, jehož součet je největší. Vydejte souřadnice začátku a konce úseku a dosažený součet.

```haskell
soucty :: Num a => [a] → (Int, Int, a)
```
Pokuste se o nějakou optimalizaci, tj. nepočítejte součty hrubou silou (zcela samostatně).

Příklad: (indexováno od 0)
```haskell
> soucty [-1,1,2,3,-4]
 (1,3,6)
```

2. Jaký význam má část ``Num a =>`` v definici funkce soucty ? Proč tam musí být?
3. Uveďte dvě možné konkrétní hodnoty proměnné a z typu funkce soucty.
4. Lze definovat ``Num a``taky pro uživatelské typy nebo musíme použít pouze předdefinované/vestavěné? Lze naši funkci soucty použít pro nějaký uživatelský typ na místě ``a`` ? (Proč ano/ne?) 

Řešení:
```haskell=
scan :: (b -> a -> b) -> b -> [a] -> [b]
scan _ acc []       = [acc]
scan f acc (x : xs) = acc : scan f (f acc x) xs

soucty :: (Num a, Ord a) => [a] -> (Int, Int, a)
soucty [] = (0, 0, 0)
soucty xs = foldr
  partSum
  (0, 0, 0)
  [ (from, to - 1, cumsum !! to - cumsum !! from)
  | from <- [0 .. length xs]
  , to   <- [from + 1 .. length xs]
  ]
 where
  cumsum = scan (+) 0 xs
  partSum :: Ord b => (a, a, b) -> (a, a, b) -> (a, a, b)
  partSum a@(_, _, x) b@(_, _, y) | x > y     = a
                                  | otherwise = b
```

### Rekonstrukce binárního stromu

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095&sid=fe143536d7d0b2e925781a412fafbdc8)


Binární vyhledávací strom je zadán jako seznam hodnot vrcholů v pořadí preorder. Definujte funkci

```haskell
readBt :: Ord a => [a] -> Bt a
```

která ze zadaného seznamu zrekonstruuje původní strom typu

```haskell
data Bt a = Void
          | Node (Bt a) a (Bt a)
```

Připomeňme, že v binárním vyhledávacím stromu platí pro každý vrchol ``v``, že všechny hodnoty v levém, resp. pravém podstromu v jsou menší, resp. větší nežli hodnota ``v``. Odtud plyne, že původní strom je zadaným seznamem určen jednoznačně.

Příklad:
```haskell
> readBt [5, 2, 4, 9]
Node (Node Void 2 (Node Void 4 Void)) 5 (Node Void 9 Void)
```

1. Definujte funkci ``readBt``.
2. Je ve vašem řešení použita nějaká funkce vyššího řádu (funkce s funkcionálními argumenty)? Pokud ne, dala by se zde nějaká smysluplně použít?
3. Je ve vašem řešení použita notace stručných seznamů (list comprehension), tj. ``[... | ...]`` ? Pokud ne, dala by se zde smysluplně použít?

Řešení:
```haskell=
data Bt a = Void
          | Node (Bt a) a (Bt a)
          deriving (Eq, Show)

readBt :: Ord a => [a] -> Bt a
readBt []       = Void
readBt (x : xs) = Node leftAns x rightAns
 where
  left     = takeWhile (<= x) xs
  right    = dropWhile (<= x) xs
  leftAns  = readBt left
  rightAns = readBt right
```

### Rostoucí posloupnosti

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089&sid=fe143536d7d0b2e925781a412fafbdc8)

Cílem je definovat funkci ``ascending,`` která na vstupu obdrží seznam hodnot (libovolného typu) a vrátí zpět seznam posloupností, který splňuje:

 - každá posloupnost je striktně rostoucí a nelze ji zleva ani zprava prodloužit
 - sloučením všech posloupností dostaneme vstupní seznam

Příklad:
```haskell
ghci> ascending [1,2,3,4,3,2,1,2]
[[1,2,3,4],[3],[2],[1,2]]

ghci> let x = [1,2,3,1,2,3] in concat (ascending x) == x
True
```

 1. Definujte typovou signaturu funkce ``ascending``.
 2. Definujte vlastní funkci.
 3. Jak byste zobecnili tuto funkci tak, aby ji bylo možné použít s libovolným porovnávacím operátorem?
 4. Bude vaše definice fungovat i na nekonečných seznamech? Pokud ano, vysvětlete proč. Pokud ne, dala by se vaše definice takto upravit? Zdůvodněte proč.

Řešení:
```haskell=
ascending :: Ord a => [a] -> [[a]]
ascending [] = []
ascending xs = ys : ascending zs where (ys, zs) = takeAscending (<) xs

takeAscending :: Ord a => (a -> a -> Bool) -> [a] -> ([a], [a])
takeAscending f []  = ([], [])
takeAscending f [x] = ([x], [])
takeAscending f (x1 : x2 : xs) | f x1 x2   = (x1 : asc, rest)
                               | otherwise = ([x1], x2 : xs)
  where (asc, rest) = takeAscending f (x2 : xs)

```

### Stromové operace

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089&sid=fe143536d7d0b2e925781a412fafbdc8)

 1. Definujte datový typ pro binární stromy.
    - Hodnoty jsou uloženy ve vnitřních uzlech.
    - Pokuste se o co nejobecnější definici.
    - Nezapomeňte na reprezentaci prázdného stromu.
 2. Definujte funkci ``replicateT``. Výsledkem ``replicateT n a`` je binární strom, který obsahuje ``n`` kopií hodnoty ``a``.
    - Výsledný strom by měl mít minimální možnou hloubku. Např. strom ``replicateT 7 a`` by měl mít hloubku 3.
 3. Definujte funkci ``zipWithT`` jako zobecnění funkce ``zipWith``. ``zipWithT f t1 t2`` sloučí prvky stromů ``t1`` a ``t2`` na stejných pozicích pomocí funkce f.
    - Pokud nemá nějaký prvek z jednoho stromu odpovídající prvek na stejné pozici v druhém stromě, tak jej do výsledného stromu nepřidávejte. Např. pro prázdný strom empty by mělo platit ``zipWithT f t empty == empty`` a ``zipWithT f empty t == empty``.
 4. Pomocí ``replicateT`` a ``zipWithT`` definujte funkci ``cut``. Funkce ``cut n t`` odstraní ze stromu ``t ``všechny vrcholy, jejichž hloubka je ostře větší než ``n``.

Řešení:
```haskell=
data Tree a = Null
            | Tree (Tree a) a (Tree a)
            deriving (Eq, Show)

replicateT :: Int -> a -> Tree a
replicateT 0 _   = Null
replicateT n val = Tree leftAns val rightAns
 where
  leftN    = (n - 1) `div` 2
  rightN   = n - 1 - leftN
  leftAns  = replicateT leftN val
  rightAns = replicateT rightN val

zipWithT :: (a -> b -> c) -> Tree a -> Tree b -> Tree c
zipWithT _ Null                  _                     = Null
zipWithT _ _                     Null                  = Null
zipWithT f (Tree leftA a rightA) (Tree leftB b rightB) = Tree leftAns
                                                              (f a b)
                                                              rightAns
 where
  leftAns  = zipWithT f leftA leftB
  rightAns = zipWithT f rightA rightB

cut :: Int -> Tree a -> Tree a
cut n tree = zipWithT const tree (replicateT ((2 ^ n) - 1) undefined)

```

### Klouzavé průměry

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078&sid=fe143536d7d0b2e925781a412fafbdc8)

Cílem je definovat binární funkci klouzave, která

- obdrží na vstupu posloupnost čísel a přirozené číslo ``n``
- a vrátí posloupnost klouzavých průměrů řádu ``n``, tj. aritmetických průměrů ``n`` sousedních prvků.

Příklad:
```haskell
klouzave [1.5, 2.5, 3.5, 4.5, 5.5] 3
[2.5,3.5,4.5]
```

1. Definujte typovou signaturu funkce ``klouzave``
2. Definujte vlastní funkci s explicitním využitím rekurze
3. Sestavte alternativní definici, tentokráte bez explicitního použití rekurze, přitom můžete využívat libovolné knihovní funkce z přiloženého seznamu.
4. Vyhýbá se alespoň jedna z vašich definic opakovaným výpočtům? Pokud ne, dala by se takto upravit? Zdůvodněte.
5. Bude některá z vašich definic fungovat i na nekonečných seznamech? Pokud ano, vysvětlete proč. Pokud ne, dala by se některá z vašich definic takto upravit? Zdůvodněte.

```haskell
take 5 $ klouzave [1..] 10
[5.5,6.5,7.5,8.5,9.5]
```

Řešení:
```haskell=
klouzave :: [Double] -> Int -> [Double]
klouzave _  0 = []
klouzave [] _ = []
klouzave xs@(first : _) n | length window < n = []
                          | otherwise         = klouzave' window mean n rest
 where
  window = take n xs
  mean   = sum window / fromIntegral n
  rest   = drop n xs

klouzave' :: [Double] -> Double -> Int -> [Double] -> [Double]
klouzave' _        mean _ []       = [mean]
klouzave' (w : ws) mean n (x : xs) = mean : klouzave' (ws ++ [x]) newMean n xs
  where newMean = (mean * fromIntegral n - w + x) / fromIntegral n


klouzave2 :: [Double] -> Int -> [Double]
klouzave2 xs n =
  [ sum window / fromIntegral n
  | d <- [0 .. length xs - n]
  , let tail   = drop d xs
  , let window = take n tail
  , length window >= n
  ]
```

### Stromovy fold

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078&sid=fe143536d7d0b2e925781a412fafbdc8)

Cílem toho problému je zobecnit funkce ``foldr`` / ``foldl`` na obecné kořenové stromy.

1. Definujte datový typ pro reprezentaci obecných kořenových stromů s ohodnocenými vrcholy:
    - snažte se o co nejobecnější definici
    - nezapomeňte na reprezentaci prázdného stromu
2. Funkce ``foldl`` a ``foldr`` zobecněte na funkci ``foldT``, která bude - namísto seznamu - procházet stromem ve vaší reprezentaci popsané v **1.**.
3. Pomocí funkce fold definujte funkci ``arita,`` která vrátí ``aritu`` (tj. maximální počet dětí přes všechny vrcholy) zadaného kořenového stromu.
4. Pomocí funkce ``foldT`` definujte funkci ``pdc``, která vrátí průměrnou délku cesty z kořene do listu (tj. součet délek všech cest z kořene do listu / počet listů).
    
Řešení:
```haskell=
-- NOTE: Tree a [] is invalid
data Tree a = Null
            | Tree a [Tree a]
  deriving (Eq, Show)

foldT :: (a -> [b] -> b) -> b -> Tree a -> b
foldT _ acc Null          = acc
foldT f acc (Tree val ts) = f val bs where bs = map (foldT f acc) ts

arita :: Tree a -> Int
arita = foldT (\_ bs -> max (length bs) (maximum bs)) 0

pdc :: Tree a -> Double
pdc tree = sumOfLengths / leafsN
 where
  (sumOfLengths, leafsN) = foldT step (0, 1) tree
  step _ bs = (sum $ map ((+ 1) . fst) bs, sum $ map snd bs)

testTree = Tree
  1
  [ Tree 1 [Null]
  , Null
  , Tree 2 [Tree 4 [Tree 5 [Null], Tree 6 [Null]]]
  , Tree 3 [Null, Null, Null, Null]
  ]
```
    
### Deleni stromu

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066&sid=fe143536d7d0b2e925781a412fafbdc7)

Rozdělte zadaný binární vyhledávací strom ``T`` na ``n+1`` binárních vyhledávacích stromů ``T_0, .. , T_n`` podle zadaných vstupních hodnot ``k_i``, ``1 <= i <= n`` tak, že ve stromě ``T_i`` jsou hodnoty ``x``, ``k_i <= x < k_i+1``, pokud jsou nerovnosti aplikovatelné.

Obrázek
![](https://i.imgur.com/sNvonyI.png)

Snažte se o efektivitu, celé podstromy patřící do jednoho pruhu zpracujte najednou.

1. Definujte datový typ pro reprezentaci binárních vyhledávacích stromů. Snažte se o co nejobecnější definici.
2. Definujte typovou signaturu funkce ``pruhy``, včetně typových tříd.
3. Funkci ``pruhy`` definujte. Budete-li používat pomocné funkce, u každé popište její význam.
4. Pokuste se stručně zdůvodnit korektnost vaší defnice.

Řešení:
[TODO]

### Run-length encoding/decoding

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066&sid=fe143536d7d0b2e925781a412fafbdc7)

Definujte funkce ``rle`` a ``rld``, které realizují run-length encoding a decoding. Funkce

```haskell
rle :: Eq a => [a] -> [Either a (a,Int)]
```

zakóduje co nejdelší úseky stejných prvků ve vstupním seznamu do dvojice ``(prvek, počet)`` typu ``Either`` s datovým konstruktorem ``Right``.

Pokud je prvek v úseku sám, kóduje se pouze prvek vnořený do typu ``Either`` s datovým konstruktorem ``Left``.

Příklad:
```haskell
> rle ”abbcccda” 
[Left 'a', Right ('b',2), Right ('c',3), Left 'd', Left 'a']
```

1. Definujte funkci ``rle`` s využitím rekurze, ale bez použití stručných seznamů či funkcí vyšších řádů (funkce s funkcionálními parametry).
2. Definujte funkci ``rle`` bez explicitního využití rekurze, ale za použití stručných seznamů či funkcí vyšších řádů.
3. Definujte typovou signaturu funkce ``rld``, která realizuje dekompresi, tj. převod ze seznamu úseků na původní seznam prvků.
4. Definujte funkci ``rld.`` Použijte přitom funkci ``map`` či ``concat``.
5. Bude některá z funkcí fungovat i na nekonečných seznamech? Proč ano nebo proč ne?

Řešení:
```haskell=
rle :: Eq a => [a] -> [Either a (a, Int)]
rle []       = []
rle (x : xs) = encoded : rle rest
 where
  ((_, count), rest) = munch x xs
  encoded            = if count == 1 then Left x else Right (x, count)

  munch x [] = ((x, 1), [])
  munch x (y : ys)
    | x /= y    = ((x, 1), y : ys)
    | otherwise = let ((_, count), rest) = munch y ys in ((x, count + 1), rest)


groups :: Eq a => [a] -> [[a]]
groups []           = []
groups xs@(x : xs') = takeWhile (== x) xs : groups (dropWhile (== x) xs')

rle2 :: Eq a => [a] -> [Either a (a, Int)]
rle2 xs =
  [ if l == 1 then Left (head group) else Right (head group, l)
  | group <- groups xs
  , let l = length group
  ]
```

### Převody mezi číselnými soustavami

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977&sid=fe143536d7d0b2e925781a412fafbdc8)

Definujte funkce:
```
prevod1 cislo puvodni
```
pro převod čísla z číselné soustavy o základu ``puvodni`` do ``dekadické číselné soustavy``, a

```
prevod2 cislo nova 
```
pro převod čísla z dekadické do číselné soustavy o základu ``nova``.
Příklad:
```haskell
> prevod1 [1,1,1,0] 2    -- převede binární 1110 do desítkové soustavy
 14
> prevod2 33 16    -- převede dekadické číslo 33 do hexadecimální soustavy
[2,1]
```

1. Doplňte typové signatury definovaných funkcí

```
prevod1 :: 
prevod2 :: 
```

2. Definujte funkci ``prevod1`` s využitím rekurze.
3. Sestavte alternativní definici funkce prevod1 s využitím alespoň jedné z funkcí ``map``, ``filter``, ``foldr`` či ``foldl``, ale bez (explicitního) použití rekurze.
4. Definujte funkci ``prevod2`` s využitím funkce ``unfold`` definované následovně:

```haskell
unfold :: (t -> Bool) -> (t -> (a, t)) -> t -> [a]
unfold done step x =  if done x then []
                                else let (y,ys) = step x 
                                     in y: unfold done step ys
```

Řešení:
```haskell=
prevod1 :: [Int] -> Int -> Int
prevod1 ds base = go ds 0
 where
  go []       acc = acc
  go (d : ds) acc = go ds (base * acc + d)

prevod1' :: [Int] -> Int -> Int
prevod1' ds base = foldl (\acc d -> base * acc + d) 0 ds


unfold :: (t -> Bool) -> (t -> (a, t)) -> t -> [a]
unfold done step x =
  if done x then [] else let (y, ys) = step x in y : unfold done step ys

moddiv a b = (a `mod` b, a `div` b)

prevod2 :: Int -> Int -> [Int]
prevod2 n base = reverse $ unfold (== 0) (`moddiv` base) n
```

### Řády prvků grupy

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977&sid=fe143536d7d0b2e925781a412fafbdc8)

Definujte unární funkci rady, která obdrží multiplikativní tabulku grupy jako matici prvků. První řádek matice obsahuje násobení grupovou jednotkou ``e`` a pořadí prvků odpovídající řádkům a sloupcům je stejné. Vydá seznam všech prvků spolu s jejich řády.

Řád prvku ``p`` je nejmenší přirozené číslo ``n ``takové, že ``n``-tá mocnina ``p`` je rovna ``e``.
1. Definujte typovou signaturu funkce ``rady``.
2. Funkci ``rady`` definujte.

Příklad:
```haskell
> rady [["e","a","b"], ["a","b","e"], ["b","e","a"]]
[("e",1), ("a",3), ("b",3)]
```
Řešení:
```haskell=
rady :: Eq a => [[a]] -> [(a, Int)]
rady []                           = []
rady table@(firstRow@(e : _) : _) = zip firstRow (map rad firstRow)
 where

  transitions =
    [ ((a, b), c)
    | (a, (b, c)) <- concat
      $ zipWith (zip . repeat) firstRow (map (zip firstRow) table)
    ]

  lookup :: Eq a => [(a, b)] -> a -> b
  lookup [] _ = error "Lookup error"
  lookup ((key, value) : ps) k | key == k  = value
                               | otherwise = lookup ps k

  mult a b = lookup transitions (a, b)
  rad x = length (takeWhile (/= e) $ iterate (mult x) x) + 1

```

### Kumulativní součty

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969&sid=fe143536d7d0b2e925781a412fafbdc8)

Je dána číselná matice ``A``. Definujte funkci

```haskell
kumulace :: Num a => [[a]] -> [[a]]
```

která z matice ``A`` vyrobí matici ``B ``stejných rozměrů (viz příklad níže).

Každý prvek na souřadnicích ``(i,j)`` bude roven součtu všech hodnot v submatici s levým horním rohem ``(0,0)``a pravým dolním rohem ``(i,j)``.

*Poznámka: Snažte se vyhnout opakování stejných výpočtů.*

Příklad:
```haskell
> kumulace[[1,1,1],[1,2,1],[0,1,0],[1,1,-4]]
[[1,2,3],[2,5,7],[2,6,8],[3,8,6]]
```

Řešení:
```haskell=
kumulace :: Num a => [[a]] -> [[a]]
kumulace []    = []
kumulace [[]]  = [[]]
kumulace table = memo
 where
  indices =
    [ [ (i, j) | (j, _) <- zip [0 ..] row ] | (i, row) <- zip [0 ..] table ]
  memo = map (map go) indices

  go (0, 0) = table !! 0 !! 0
  go (i, 0) = memo !! (i - 1) !! 0 + table !! i !! 0
  go (0, j) = memo !! 0 !! (j - 1) + table !! 0 !! j
  go (i, j) =
    (memo !! (i - 1) !! j)
      + (memo !! i !! (j - 1))
      - (memo !! (i - 1) !! (j - 1))
      + (table !! i !! j)
```

### Doplnění hypergrafu

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969&sid=fe143536d7d0b2e925781a412fafbdc8)

Hypergraf je zadán množinou vrcholů a množinou hyperhran, což jsou alespoň dvouprvkové podmnožiny množiny vrcholů. Naší cílem je definovat funkci doplnění, která doplní do hypergrafu ``H`` všechny dvouprvkové (hyper)hrany pro ty dvojice vrcholů, které nejsou společně obsaženy v žádné hyperhraně vstupního hypergrafu ``H``. Funkce tedy např. z hypergrafu s vrcholy ``{1,2,3,4,5}`` a hyperhranani ``{1,3,5}`` a ``{2,3,4}`` vytvoří hypergraf se stejnými vrcholy a hyperhranami ``{1,3,5},{2,3,4},{1,2},{1,4},{5,2}`` a ``{5,4}``

1. Definujte datový typ pro reprezentaci hypergrafu. Pokuste se o co nejobecnější definici (vrcholy mohou být reprezentovány nejen čísly, ale i znaky, řetězci apod.)
2. Specifikujte typovou signaturu funkce

```
doplneni ::
```

3. Funkci definujte.
    
Řešení:
```haskell=
data HGraph a = HGraph [a] [[a]]
  deriving (Eq, Show)

fromMaybe :: a -> Maybe a -> a
fromMaybe def Nothing  = def
fromMaybe _   (Just a) = a

lookup' :: Eq a => [(a, b)] -> a -> Maybe b
lookup' [] _ = Nothing
lookup' ((key, value) : ps) k | key == k  = Just value
                              | otherwise = lookup' ps k

doplneni :: Eq a => HGraph a -> HGraph a
doplneni (HGraph vs es) = HGraph vs (es ++ newEdges)
 where
  trans_table = zip [0 ..] vs
  newEdges =
    [ [ fromMaybe undefined (lookup' trans_table i)
      , fromMaybe undefined (lookup' trans_table j)
      ]
    | i <- [0 .. length vs - 1]
    , j <- [0 .. i - 1]
    , not $ isEdge i j
    ]

  isEdge i j = or
    [ elem (fromMaybe undefined (lookup' trans_table i)) edge
        && elem (fromMaybe undefined (lookup' trans_table j)) edge
    | edge <- es
    ]

```

### Analýza textu (5 bodů)

Zdroj: [Zkouška 10.6.2019 (Dvořák + Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11954)

Na vstupu je zadán text jako hodnota typu String. Naším cílem je definovat binární funkci stat text n, která

    obdrží takový text a přirozené číslo n
    vrátí všechna slova z tohoto textu o délce alespoň n, setříděná lexikograficky
    každé slovo s čísly řádků, kde se slovo vyskytuje

Řádky jsou ukončeny znakem '\n'. Slovo je každý maximální podřetězec textu neobsahující mezeru ' ', tabulátor '\t' či konec řádku '\n'.
(a) Definujte datovou strukturu pro reprezentaci oboru hodnot funkce stat (pomocí data nebo type).
(b) Definujte typovou signaturu funkce stat s použití datové struktury z (a).
(c) Funkci stat definujte.

### Označkování stromu

Zdroj: [MFF Forum: Zkouška 21.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11747)

Ohodnotit vrcholy obecného n-arní stromu v pořadí v jakém jsou vrcholy uzavírány, takže post-fixově. Měla se napsat datová struktura pro strom, typová hlavička fce a pak tu funkci implementovat:

Kód: Vybrat vše

data Strom a = Nil | Strom a [Strom a]
label::Strom a -> Strom (a,Int)


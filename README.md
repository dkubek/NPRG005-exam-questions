# NeprocProg Cure

[TOC]

## Prolog

### Topologické uspořádání grafu

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095)

Je dán orientovaný graf G pomocí seznamů sousedů. Zjistěte, jestli lze graf G
topologicky uspořádat a pokud ano, vydejte seznam vrcholů v topologickém
pořadí.

Příklad:

```prolog
?- topo([a-[],b-[a,c],c-[a],d-[a,c]],Usp).
Usp = [b,d,c,a]
```

 1) Definujte příslušný predikát topo/2 v jazyce Prolog.
 2) Odhadněte časovou složitost vašeho řešení. Odhad zdůvodněte.
 3) Jsou některé z vašich predikátů koncově rekurzivní ? Pokud ano, vysvětlete,
 které to jsou, a jaký to má význam. Pokud ne, vysvětlete, zdali by se dal
 některý takto upravit.

Řešení:

```prolog
remove_vertex([], _, []).
remove_vertex([ Vertex-_ | Graph ], Vertex, Out) :-
    remove_vertex(Graph, Vertex, Out),
    !.
remove_vertex([ V-Ns | Graph ], Vertex, [ V-NewNs | Ans ]) :-
    remove_vertex(Graph, Vertex, Ans),
    (
        member(Vertex, Ns)
    ->
        select(Vertex, Ns, NewNs)
    ;
        NewNs = Ns
    ).

topo(Graph, Usp) :-
    topo_(Graph, [], Usp).

topo_([], Acc, Acc).
topo_(Graph, Acc, Out) :-
    member(Min-[], Graph),
    remove_vertex(Graph, Min, NewGrap),
    topo_(NewGrap, [Min | Acc], Out).
```

### Diskrepanční vrstvy

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095)

Napište predikát ``diskr/2``, který dostane binární strom (s konstruktory
``t/3`` a ``nil/0``) a vrátí seznam seznamů vrcholů stromu, kde v jednom
vnitřním seznamu jsou všechny vrcholy, ke kterým se při průchodu od kořene
dostaneme se stejným počtem kroků doprava. Vnější seznam je od nejlevější
vrstvy, na pořadí ve vnitřních seznamech nezáleží.

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
2. Je ve vašem řešení použit řez (!) nebo negace? Pokud ano, změní se něco,
   když řez / negaci vypustíme? Pokud ne, dal by se řez / negace někde
   smysluplně využít?
3. Lze u predikátu ``diskr/2`` obrátit směr výpočtu? Podrobněji: dle příkladu
   předpokládáme volání diskr(+,-). Bude fungovat i volání diskr(-, +), tj.
   zadáme seznam diskrepančních vrstev, a na výstupu obdržíme strom?
   Vysvětlete.

Řešení:

```prolog
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

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089)

Cílem úlohy je definovat predikát allTrees, který pro daný seznam hladin
vygeneruje všechny možné binární stromy.

- Hladinou rozumíme seznam prvků, které se nacházejí ve stejné hloubce
- Můžete předpokládat, že každá hladina má nanejvýš dvojnásobek prvků předchozí
  hladiny (ale může jich mít méně).
- Hladiny vygenerovaného stromu musejí odpovídat hladinám specifikovaných ve
  vstupním seznamu.

Např. pro seznam ``[[1],[2,3],[4]]`` dostaneme následující 4 stromy:

```none
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
 4. Lze vaší definici použít opačným směrem? Tj. nalezne váš predikát seznam
    hladin pokud specifikujete pouze výsledný strom? Vysvětlete.

Řešení:

```prolog
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

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089)

Je zadán neorientovaný graf *G* a množina vrcholů *M*. Zjistěte, zda *M* a
doplněk *M* tvoří bipartitní rozklad grafu *G* (tj. každá hrana grafu má právě
jeden koncový vrchol v množině *M*). Pokud ano, vydejte druhou množinu
rozkladu.

```prolog
?- bip([a-[c,d], b-[d], c-[a], d-[a,b]], [a,b], D).
    D = [c,d]

?- bip([a-[c,d], b-[d], c-[a], d-[a,b]], [b,c], D).
    false
```

 1. Definujte predikát ``bip/3``.
 2. Napište o jednotlivých predikátech ve vašem řešení, zda jsou koncově rekurzivní.

Řešení:

```prolog
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

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078)

Truhlář má dostatek trámů délky ``D`` a seznam ``Xs`` délek trámů, které
potřebuje nařezat. V seznamu ``Xs`` se délky mohou opakovat.

Cílem problému je sestavit predikát ``rezy(+D, +Xs, -N, -Vss)``, který

- rozdělí požadované délky do skupin, které se mají nařezat z jednoho trámu
- truhlář přitom používá hladový algoritmus, tj. pro každou délku použije první
  trám, z něhož lze ještě požadovanou délku odřezat
- vrátí celkový počet řezaných trámů N
- a seznam seznamů Vss (délky N), jehož každý prvek reprezentuje dělení jednoho
  trámu (případný zbytek se neuvádí).

```prolog
?- rezy(5,[3,2,2,2,2,1,4], N, V).
N=4, V=[[3,2],[2,2,1],[2],[4]]
```

1. Definujte predikát ``rezy/4.`` Definice případných pomocných predikátů
   prosím opatřete vysvětlujícím komentářem.
2. Je některý z vašich predikátů koncově rekurzivní? Pokud ano, vysvětlete,
   který to je a jaký to má význam.
3. Pokud ne, dal by se některý takto upravit? Odpověď prosím zdůvodněte.

Řešení:

```prolog
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

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078)

Je zadán seznam množin ``Mss``. Chceme všemi možnými způsoby vybrat a vrátit v
seznamu reprezentanty daných množin v odpovídajícím pořadí s podmínkou, že
konkrétní reprezentanti v jednom výběru jsou různí.

Příklad:

```prolog
?- reprezentanti([[1],[1,2,3],[1,3,4]], R).
R = [[1,2,3],[1,2,4],[1,3,4]]
```

1. Sestavte predikát ``reprezentanti(+Mss, -Rss)``.
2. Stručně vysvětlete, proč je vaše definice korektní.
3. Je ve vašem programu použit řez ``(!)`` ? Jde o řez červený (mění
   deklarativní význam programu) či zelený (nemění d.v.)? Pokud ne, je řez
   nezbytný pro definici některého vestavěného predikátu / operátoru, který
   jste ve vašem řešení použili? Jde o řez červený (mění deklarativní význam
   programu) či zelený (nemění d.v.)?

```prolog
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

### Hammerstein

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066)

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

1. Doplňte jednu (opravdu jen jednu) chybějící klauzuli za uvedené pravidlo
   tak, aby výsledná procedura korektně setřídila vstupní seznam přirozených
   čísel. Na výstupu bychom měli obdržet jen jediné řešení.
2. V definici pravidla je použit řez (!). Jde o zelený (nemění deklarativní
   význam) či červený řez (mění d.v.) ? Vysvětlete! Obsahuje některá z vašich
   klauzulí, (doplněná v(a) nebo (b)) zelený či červený řez?
3. Jaký známý třídící algoritmus výše uvedený kód implementuje? Pokud neznáte
   název, můžete alespoň slovně popsat, jak ``setrid/2`` funguje.
4. *VOLITELNE*: Lze u procedury ``setrid/2`` obrátit směr výpočtu?

```prolog
setrid(-Xs,+Ys) :- Xs je seznam přirozených čísel ze seznamu Ys setříděný vzestupně
```

Pokud ne, šel by kód jednoduše upravit tak, aby se výsledný predikát
(pojmenovaný třeba ``setrid2/2``) dal korektně volat oběma způsoby?

Řešení:

```prolog
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

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066)

Do země Mobilia, v níž je každý občan vybaven chytrým telefonem, přicestoval
Cestovatel, nakažený virovým onemocněním. Všichni ostatní byli přitom ještě
zdraví. Můžeme předpokládat, že virus se přenese z jedné osoby na druhou, pokud
spolu strávili ve vzdálenosti menší než 2m alespoň čas ``K``, kde ``K`` je
známá kritická hodnota. Díky chytrým telefonům máme pro každého občana Mobilie
seznam záznamů jeho kontaktů, kde každý takový záznam pro osobu ``A`` obsahuje
identifikaci osoby ``B``, která se k němu přiblížila do vzdálenosti ``< 2m``
čas setkání a délku setkání.

Cílem je sestavit program, který na základě takových záznamů vrátí seznam
infikovaných osob.

1. V jazyce Prolog popište datovou strukturu pro reprezentaci jednoho záznamu
   kontaktu občana Mobilie popsaného výše.
2. V jazyce Prolog navrhněte reprezentaci položek VstupníhoSeznamu, přičemž
   každá položka bude obsahovat indentifikaci občana Mobilie a seznam záznamů
   jeho kontaktů.
3. Sestavte predikát ``inf/4``, který obdrží

    ```none
    VstupníSeznam
    identifikaci Cestovatele
    kritickou hodnotu K
    ```

    a vrátí seznam infikovaných.

    U každého pomocného predikátu prosím v poznámce popište jeho význam.

    *Volitelné:* výstupní seznam můžete uspořádat dle délky kontaktu s infikovanými
    do nerostoucí posloupnosti.

4. Odhadněte časovou složitost vašeho řešení.
5. Je některý z vašich predikátů koncově rekurzivní ? Pokud ano, vysvětlete,
   jaký to má význam. Pokud ne , dal by se některý takto upravit?

Řešení:

```prolog
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

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977)

Definujte binární predikát ``aspon2/2``, který

- obdrží seznam výrokových proměnných (reprezentovaných atomy), v němž je každá
  proměnná ohodnocena hodnotou true nebo false
- vrátí seznam všech takových ohodnocení týchž proměnných, v němž se každé
  ohodnocení bude od vstupního lišit v hodnotách alespoň 2 proměnných.

Příklad:

```prolog
?- aspon2([x1-true, x2-false, y-true], V).
  V =  [  [x1-false, x2-true, y-true],
          [x1-false, x2-false, y-false],
          [x1-true, x2-true, y-false],
          [x1-false, x2-true, y-false] ]
```

Řešení:

```prolog
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

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977)

Graf je zadán jako seznam svých vrcholů se seznamy sousedů (viz příklad).
Definujte binární predikát ``troj(+Graf, -SeznamTrojuhelniku)`` který k
takovému grafu vrátí seznam všech jeho trojúhelníků. Ve výsledném seznamu by se
každý trojúhelník měl vyskytovat právě jednou (``t(a,b,c)``, ``t(b,c,a)`` a
``t(c,a,b)`` jsou stejné trojúhelníky).

Příklad:

```prolog
?- troj([a-[b,c,d],b-[a,c],c-[a,b,d],d-[a,c],e-[]], S).
     S = [t(a,b,c), t(a,c,d)]
```

Řešení:

```prolog
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

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969)

Formule výrokového počtu jsou sestavené z (výrokových) proměnných ve funktoru
``var/1`` a logických spojek negace, konjunkce a disjunkce (bez konstant). Dále
máte dány v argumentech predikátu ``gen/3`` číslo ``k`` pro velikost formule a
seznam jmen proměnných. Generujte backtrackingem všechny logické formule
(každou jednou), které obsahují proměnné ze seznamu a ve kterých je počet
spojek a výskytů proměnných dohromady právě ``k``.

Definujte predikát ``gen(+K, +Jmena, -Fle)``. Na pořadí generovaných formulí
nezáleží, ale měli byste vygenerovat každou právě jednou. K řešení není potřeba
predikát ``=../2`` (univ).

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

```prolog
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

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969)

Jeden objekt je zadán uspořádaným seznamem dvojic klíč-hodnota. Na vstupu máte
seznam objektů. Napište proceduru ``koncept/2``, která vyrobí nejmenší koncept
zahrnující všechny vstupní objekty. Koncept je seznam dvojic
``klíč-seznam_hodnot``. Koncept zahrnuje objekt, pokud koncept má všechny klíče
objektu a v seznamu hodnot příslušného klíče u konceptu je obsažena hodnota
klíče u objektu. Pokud objekt nějaký klíč konceptu nemá, bude v seznamu hodnot
konceptu hodnota ``nedef``.

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

```prolog
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

Máte dány dva řetězce, u kterých nevíte jejich vzájemnou orientaci. Najděte a
vydejte v seznamu všechna jejich vzájemná neprázdná překrytí.

Příklad:

```prolog
?- prekryti([a,t,c,t,c],[c,t,c,c], V).
V = [a,t,c,t,c,t,c,c],[a,t,c,t,c,c],[a,t,c,t,c,c,t,c]]
```

Řešení:

```prolog
id_or_reverse(X, X).
id_or_reverse(X, Y) :-
    reverse(X, Y),
    X \= Y.

prekryti(Xs, Ys, Out) :-
    prekryti_(Xs, Ys, [], Out).

prekryti_(Xs, Ys, Acc, Out) :-
    is_prekryti(Xs, Ys, P),
    \+ member(P, Acc),
    prekryti_(Xs, Ys, [P | Acc], Out),
    !.
prekryti_(_, _, Acc, Acc) :- !.

is_prekryti(Xs, Ys, Out) :-
    id_or_reverse(Xs, X),
    id_or_reverse(Ys, Y),

    append(_, BodyTailX, X),
    append(HeadBodyY, TailY, Y),

    BodyTailX = HeadBodyY,
    BodyTailX \= [],

    append(X, TailY, Out).
```

### Neporovnatelné prvky částečně uspořádané množiny

Zdroj: [MFF Forum: Zkouška 10.6.2019 (Dvořák + Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11954)

Částečně uspořádaná množina je popsána seznamem termů tvaru ``x -> y`` s
významem x pokrývá y (tj. ``x > y`` a současně ``x ≥ z ≥ y`` implikuje ``x =
z`` nebo ``y = z``).

Definujte predikát ``nepor/2``, který k takto zadané množině vrátí seznam všech
dvojic vzájemně neporovnatelných prvků (tj. dvojic ``x``,``y`` takových, že
neplatí ``x ≥ y`` ani ``x ≤ y``).

Příklad:

```prolog
?- nepor([a->b, a->c, b->d, e->f], N).
N = [a-e,a-f,b-c,b-e,b-f,c-d,c-e,c-f,d-e,d-f]
```

Řešení:

```prolog
ge(_, X, X).
ge(R, X, Y) :- member(X -> Y, R).
ge(R, X, Y) :- member(X -> Z, R), ge(R, Z, Y).


collect_variables(Rel, Vars) :-
    collect_variables_(Rel, [], Vars).

collect_variables_( [], Acc, Ans) :-
    sort(Acc, Ans).
collect_variables_( [X -> Y | Rs], Acc, Ans) :-
    collect_variables_(Rs, [X, Y | Acc], Ans).

pair(X, Y, X-Y).

pairs(Vars, Pairs) :-
    select(Var, Vars, RestVars),
    !,
    maplist(pair(Var), RestVars, Tmp),
    pairs(RestVars, Ans),
    append(Tmp, Ans, Pairs).
pairs([], []).

is_nepor(R, X-Y) :-
    \+ ge(R, X, Y),
    \+ ge(R, Y, X).

nepor(Rel, Out) :-
    collect_variables(Rel, Vars),
    pairs(Vars, Pairs),
    include(is_nepor(Rel), Pairs, Out).
```

### Lexikograficky předchozí permutace

Zdroj: [MFF Forum: Zkouška 21.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11747)

Nalezněte lexikograficky předchozí permutaci. Pokud neexistuje tak ``false``.

Řešení:

```prolog
?- prev([1,2,6,3,4,5,7],V).
V = [1,2,5,7,6,4,3]
```

```prolog
find_longest_ascending([], [], []).
find_longest_ascending([X], [X], []).
find_longest_ascending([X1, X2 | Xs], [X1], [X2 | Xs]) :-
    X1 > X2,
    !.
find_longest_ascending([X1, X2 | Xs], [ X1 | Ans ], Rest) :-
    X1 < X2,
    find_longest_ascending([X2 | Xs], Ans, Rest).

replace([], _, _, []).
replace([X | Xs], X, Y, [Y | Xs]) :- !.
replace([R | Xs], X, Y, [R | Ans]) :-
    replace(Xs, X, Y, Ans).

prev(Perm, Prev) :-
    find_longest_ascending(Perm, Asc, Rest),
    reverse(Asc, Rev),
    member(X, Rev),
    Y is X - 1,
    member(Y, Rest),
    !,
    replace(Asc, X, Y, NewAsc),
    replace(Rest, Y, X, NewRest),
    reverse(NewRest, FinalRest),
    append(NewAsc, FinalRest, Prev).

```

### Frekvence

Zdroj: [MFF Forum: Zkouška 26.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11756)

Definujte predikát ``frekv/2``, který obdrží seznam konstant
a vrátí frekvence dvojic za sebou jdoucích konstant.
Výstupní reprezentaci si zvolte a popište pro vzorový vstup.

```prolog
?- frekv([a,b,a,b,c], P).
P = [f(a-b,2), f(b-a,1), f(b-c,1)]
```

Řešení:

```prolog
frekv(List, Freq) :-
    frekv_(List, [], Freq),
    !.

frekv_([], Acc, Acc).
frekv_([_], Acc, Acc).
frekv_([ X1, X2 | Xs ], Freq, Ans) :-
    increase_frequency(X1-X2, Freq, NewFreq),
    frekv_([ X2 | Xs ], NewFreq, Ans).


increase_frequency(X-Y, Freq, [f(X-Y, NewN) | Rest]) :-
    select(f(X-Y, N), Freq, Rest),
    NewN is N + 1,
    !.
increase_frequency(X-Y, Freq, [f(X-Y, 1) | Freq]).
```

### Časové ohodnocení DFS

Zdroj: [MFF Forum: Zkouška 26.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11756)

Je dán orientovaný acyklický graf (DAG) o ``n`` vrcholech pomocí seznamu
sousedů. Procedura ``dfs/2`` projde graf do hloubky a přidá k vrcholům časy
otevření a uzavření v rozsahu od ``1`` do ``2n``. Na pořadí vrcholů na výstupu
nezáleží:

Definujte predikát ``dfs/2`` a napište konkrétní výstup vašeho programu na
vzorovém grafu z příkladu níže.

Příklad:

```prolog
?- dfs([c-[d], a-[b,c], b-[d,e], d-[], e-[]], V).
V = [v(a,1,10,[b,c]), v(c,2,5,[d]), v(d,3,4,[]), v(b,6,9,[e]), v(e,7,8,[])]
```

Řešení:

```prolog
dfs(Graph, Out) :-
    member(Start-_, Graph),
    Stack = [ Start ],
    Opened = [],
    Closed = [],
    Time = 1,
    dfs_(Graph, Stack, Time, Opened, Closed, Out).

dfs_(_, [ ], _, _, Out, Out).
dfs_(Graph, [ Vertex | Vs ], Time, Opened, Closed, Out)  :-
    member(t(Vertex, _, _), Closed),
    dfs_(Graph, Vs, Time, Opened, Closed, Out),
    !.
dfs_(Graph, [ Vertex | Vs ], Time, Opened, Closed, Out) :-
    NewTime is Time + 1,
    (
        select(Vertex-InTime, Opened, NewOpened)
    ->
        dfs_(
            Graph,
            Vs,
            NewTime,
            NewOpened,
            [t(Vertex, InTime, Time) | Closed],
            Out
        )
    ;
        member(Vertex-Neigbours, Graph),
        append(Neigbours, [Vertex | Vs], NewStack),

        dfs_(
            Graph,
            NewStack,
            NewTime,
            [ Vertex-Time | Opened ],
            Closed,
            Out
        )
    ).
```

### Splay

Zdroj: [MFF Forum: Zkouška 6. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10961)

Naprogramujte predikát ``splay(+Hodnota, +BinarniVyhledavaciStrom,
-Vysledek)``, který provede funkci ``splay`` (přesune daný vrchol až do kořene
pomoci rotací) na ``Hodnotu``. Pokud ``Hodnota`` ve stromě není, pak se splay
provede na bezprostredního předchůdce/následníka.

```prolog
TestTree = tree(
    tree(
        tree(
            tree(null, 1, null),
            2,
            tree(null, 3, null)
            ),
        4,
        tree(null, 5, null)
        ),
    6,
    tree(
        tree(null, 7, null),
        8,
        tree(null, 9, null)
        )
    ).
```

Řešení:

```prolog
splay(X, T, T) :-
    T = tree(_, X, _),
    !.
splay(X, T, Out) :-
    T = tree(Left, Y, Right),
    (
        X < Y
    ->
        splay(X, Left, Ans),
        Ans = tree(LeftAns, Z, RightAns),
        Out = tree(LeftAns, Z, tree(RightAns, Y, Right))
    ;
        splay(X, Right, Ans),
        Ans = tree(LeftAns, Z, RightAns),
        Out = tree(tree(Left, Y, LeftAns), Z, RightAns)
    ),
    !.
splay(_, T, T) :- T = tree(null, _, null).
```

### Skládání konstantních úseků

Zdroj: [MFF Forum: Zkouška 6. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10961)

Na vstupu máme seznam po částech konstantních funkcí ``Fs``, kde každá funkce
je ve tvaru ``DelkaUseku-Hodnota``. Všechny funkce začinají v ``0`` a po konci
posledního useku pokračují hodnotou ``0``. Máme vytvořit nejmenší novou funkci
takovu, že v každém bodě je větší rovna všem zadaným funkcím.

Příklad:

Dvě funkce: první má na intervalu ``[0, 2)`` hodnotu 5, na intervalu ``[2, 4)``
hodnotu 3 a na intervalu ``[4, inf)`` hodnotu 0. Druhá má na intervalu ``[0,
3)`` hodnotu 4 a na intervalu ``[3, inf)`` hodnotu 0.
Vysledkem je funkce ``[2-5, 1-4, 1-3]``.

```prolog
?- combine([[2-5, 2-3], [3-4]], G)
G = [2-5, 1-4, 1-3]
```

Řešení:

```prolog
combine([], []).
combine([ Base | Fs ], G) :-
    combine_(Fs, Base, G).

combine_([], Base, Base).
combine_([F | Fs], Base, Out) :-
    extend_base(Base, F, NewBase),
    combine_(Fs, NewBase, Out),
    !.

extend_base(Base, [], Out) :-
    merge_adjecent(Base, Out).
extend_base([], F, F).
extend_base([ Length-BaseVal | BLVs ], [ Length-FVal | FLVs], Out) :-
    extend_base(BLVs, FLVs, Ans),
    (
        BaseVal > FVal
    ->
        Out = [Length-BaseVal | Ans]
    ;
        Out = [Length-FVal | Ans]
    ).
extend_base([ BaseLength-BaseVal | BLVs ], [ FLength-FVal | FLVs], Out) :-
    (
        BaseLength < FLength
    ->
        Remainder is FLength - BaseLength,
        extend_base(
            [ BaseLength-BaseVal | BLVs ],
            [ BaseLength-FVal, Remainder-FVal | FLVs],
            Out
        )
    ;
        Remainder is BaseLength - FLength,
        extend_base(
            [ FLength-BaseVal, Remainder-BaseVal | BLVs ],
            [ FLength-FVal | FLVs],
            Out
        )
    ).

merge_adjecent([], []).
merge_adjecent([P], [P]).
merge_adjecent([Length1-Val, Length2-Val | LVs], Out) :-
    NewLength is Length1 + Length2,
    merge_adjecent([NewLength-Val | LVs], Out),
    !.
merge_adjecent([X, Y | LVs], [X | Out]) :-
    merge_adjecent([Y | LVs], Out).
```

### Kružnice v grafu

Zdroj: [MFF Forum: Zkouška 22.6.](http://forum.matfyz.info/viewtopic.php?f=169&t=11412)

Máme daný orientovaný graf reprezentovaný jako ``[vrchol-[seznam
sousedů]|...]``, zjistěte, zda v něm je orientovaná kružnice, a pokud ano,
vraťte vrcholy nějaké takové kružnice v tom pořadí, jak jsou na kružnici. Chce
se polynomiální řešení.

Příklad:

```prolog
?- cycle([a-[b,c,d],b-[c],c-[a,b,d],d-[a,c],e-[]], C)
C = [a, c, b]
```

Řešení:

```prolog
cycle(Graph, Cycle) :-
    member(Start-_, Graph),
    Stack = [Start],
    Path = [],
    cycle_(Graph, Stack, Path, Cycle),
    !.

cycle_(Graph, [ Vertex | Vs ], Path, Out) :-
    (
        member(Vertex, Path)
    ->
        append(Cs, [Vertex | _], Path),
        Out = [Vertex | Cs]
    ;
        member(Vertex-Neighbours, Graph),
        append(Neighbours, Vs, NewStack),
        cycle_(Graph, NewStack, [Vertex | Path], Out)
    ).
```

### Vypustění nejvýše dvou prvků

Zdroj: [MFF Forum: Zkouska 20.9.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11466)

Definujte predikát ``dif/2``, který obdrží seznam ``S``, a vrátí seznam všech
seznamů, které obdržíte z ``S`` vypuštěním nejvýše 2 prvků. Pořadí prvků ve
výstupních seznamech se nemění.

Příklad:

```prolog
?- dif([a,b,c],D).
D=[[a,b,c], [b,c], [a,c], [a,b], [a], [c]]
```

Řešení:

```prolog
smaller_than(N, Xs) :- length(Xs, K), K < N.

dif(List, Out) :-
    length(List, N),
    MinSize is N - 2,
    dif_(List, Ans),
    exclude(smaller_than(MinSize), Ans, Out).

cons(X, Xs, [X | Xs]).

dif_([], [[]]).
dif_([X | Xs], Out) :-
    dif_(Xs, Tmp),
    maplist(cons(X), Tmp, Appended),
    append(Appended, Tmp, Out).
```

### Vrcholové pokrytí minimální k inkluzi

Zdroj: [MFF Forum: Zkouska 20.9.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11466)

Je zadán neorientovaný graf ``G`` a jeden jeho vrchol ``v``. Definujte predikát
``pokryti/3,`` který postupně vrátí všechna vrcholová pokrytí grafu ``G``,
která obsahují vrchol ``v`` a jsou minimální vzhledem k inkluzi.

Množina vrcholů ``V`` grafu je vrcholovým pokrytím, pokud každá hrana má
alespoň jeden vrchol v množině ``V``.

1. Na příkladě popište, jakou reprezentaci grafu budete používat.
2. Definujte predikát
```pokryti(+Graf, +Vrchol, -VPokrytí)```
kde Graf je zadán v reprezentaci popsané v **1.)**.

Řešení:

```prolog
collect_nodes(Graph, Nodes) :-
    collect_nodes(Graph, [], NodesDup),
    sort(NodesDup, Nodes),
    !.

collect_nodes([], Acc, Acc).
collect_nodes([Node-Neighbours | Ns], Acc, Ans) :-
    append([Node | Neighbours], Acc, NewAcc),
    collect_nodes(Ns, NewAcc, Ans).


covers(Cover, U, V) :-
    member(U, Cover);
    member(V, Cover).

is_cover([], _).
is_cover([ Vertex-Neighbours | Rest ], Cover) :-
    maplist(covers(Cover, Vertex), Neighbours),
    is_cover(Rest, Cover).


pokryti(Graph, Vertex, Cover) :-
    collect_nodes(Graph, Nodes),
    select(Vertex, Nodes, Rest),
    !,
    pokryti_(Graph, Vertex, Rest, Cover).

pokryti_(Graph, Vertex, [], [Vertex]) :-
    is_cover(Graph, [Vertex]).
pokryti_(Graph, Vertex, Nodes, Cover) :-
    (
        select(_, Nodes, Rest),
        is_cover(Graph, [Vertex | Rest])
    ->
        pokryti_(Graph, Vertex, Rest, Cover)
    ;
        Cover = [Vertex | Nodes]
    ).
```

### !! Rozděl

Zdroj: [MFF Forum: Zkouška 13. 9. 2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11457)

Napište predikát ``rozdel(+Mnozina,-Rozdeleni)``, který rozdělí množinu na
neprázdné podmnožiny. Všechny možnosti rozdělení pak vrátí spojené v jednom
seznamu.

Příklad:

```prolog
?- rozdel([a,b,c],X).
X = [[a, b, c], [[a, b], [c]], [[a], [b, c]], [[a, c], [b]], [[a], [b], [c]]].
```

???

### Nezávislé množiny

Zdroj: [MFF Forum: Zkouška 13. 9. 2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11457)

Napište predikát ``nez(+Graf,+Vrchol.-NezMn)``, který vezme graf a jeden jeho
zadaný vrchol a postupně vydává všechny jeho největší nezávislé množiny
obsahující daný vrchol.

Příklad:

```prolog
nez(g([a,b,c,d,e],[a-b,b-c,b-d,c-d]),a,X).
X=[a,c,e];
X=[a,d,e].
```

Řešení:

```prolog
is_edge(g(_, Edges), U, V) :-
    member(U-V, Edges), !;
    member(V-U, Edges), !.

is_not_edge(Graph, U, V) :- \+ is_edge(Graph, U, V).

is_independent(_, []) :- !.
is_independent(Graph, [ V | Vs ]) :-
    maplist(is_not_edge(Graph, V), Vs),
    is_independent(Graph, Vs).

is_subset([], []).
is_subset([ X | Xs ], [ X | Ys ]) :-
    is_subset(Xs, Ys).
is_subset([_ | Xs], Ys) :-
    is_subset(Xs, Ys).

nez(Graph, Vertex, MaxIndSet) :-
    g(Vertices, _) = Graph,
    length(Vertices, N),
    nez_(Graph, N, MaxIndSet),
    !,
    member(Vertex, MaxIndSet).

nez_(Graph, N, IndSet) :-
    N > 0,
    Graph = g(Vertices, _),

    length(IndSet, N),
    is_subset(Vertices, IndSet),
    is_independent(Graph, IndSet).
nez_(Graph, N, IndSet) :-
    NewN is N - 1,
    NewN > 0,
    nez_(Graph, NewN, IndSet).
```

### Cykly délky alespoň N

Zdroj: [MFF Forum: Zkouška 6. 6. 2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11380)

Na vstupu máme graf reprezentovaný jako

```prolog
graf(SeznamVrcholu, SeznamHran)
```

(bylo ale dovoleno si reprezentaci grafu změnit) a číslo ``N``. Máme určit,
jestli v grafu existuje cyklus délky alespoň ``N``. Pokud ano, program alespoň
jeden takový cyklus vypíše, pokud ne, vrátí fail.

*Pozn.: Problém je NP-úplný, tzn. očekává se řešení typu hrubá síla.*

Řešení:

```prolog
subsets([], []).
subsets([ H | T ], [ H | Out ]) :-
    subsets(T, Out).
subsets([ _ | T], Out) :-
    subsets(T, Out).

is_edge(graph(_, Edges), U, V) :-
    member(U-V, Edges), !;
    member(V-U, Edges), !.

is_path(Graph, [], Start, End) :-
    is_edge(Graph, Start, End).
is_path(Graph, Vertices, Start, End) :-
    select(Vertex, Vertices, Rest),
    is_edge(Graph, Start, Vertex),
    is_path(Graph, Rest, Vertex, End).

is_cycle(Graph, Vertices) :-
    select(Start, Vertices, Rest),
    is_path(Graph, Rest, Start, Start),
    !.

cycle_n(Graph, N, Cycle) :-
    graph(Vertices, _) = Graph,
    length(Vertices, K),
    K >= N,
    !,

    between(N, K, N_),
    length(Cycle, N_),
    subsets(Vertices, Cycle),
    is_cycle(Graph, Cycle),
    !.
```

### Termy

Zdroj: [MFF Forum: Zkouška 29.5.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11357)

Sestavte predikát ``termy/1``, který postupně vrací termy složené z funktorů
``bin/2``, ``un/1`` a ``const/0``. Výstupem bude tedy korektně sestavený term.
Predikát by měl postupně vrátit všechna řešení, sice v libovolném pořadí, ovšem
každé právě jednou.

Příklad:

```prolog
?- termy(V).
V=const;
V=un(const);
V=bin(const,const);
V=un(un(const));
V=un(bin(const,const));
V=bin(un(const),un(const));
```

Řešení:

```prolog
termy(V) :-
    length(Slots, _),
    termy_(Slots, V).

termy_([_], const).
termy_([_ | Slots], un(T)) :-
    termy_(Slots, T).
termy_([_ | Slots], bin(T1, T2)) :-
    append(S1, S2, Slots),
    termy_(S1, T1),
    termy_(S2, T2).
```

### !! Porovnání multimnožin

Zdroj: [MFF Forum: Zkouška 29.5.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11357)

Multimnožinu lze specifikovat seznamem termů ``Prvek-Pocet``. Sestavte predikát
``mensi/2``, který porovná multimnožiny ``A`` a ``B`` následovně:

- ``mensi(A,B)`` je ``true`` právě tehdy, pokud v ``B`` existuje nějaký prvek,
  co není v ``A`` takový, že je větší než všechny prvky z ``A``, které nejsou
  v ``B``.

```prolog
?- mensi([c-3,b-2,a-1],[d-1,b-3])
true
?- mensi([c-3,b-2,a-1],[c-1,b-3])
fail
```

Řešení:
???

### Plánování výroby

Zdroj: [MFF Forum: Zkouška 13. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10963)

Napište predikát, který naplánuje pokud možno optimální (nutné použít nějakou
jednoduchou heuristiku) rozvrh výroby na strojích. Na vstupu je seznam délek
operací (např. ``[3,3,2,6,4]``) a maximální čas běhu (např. ``10``). Operace je
možné plánovat na paralelně běžící stroje, chceme, aby celkový počet potřebných
strojů byl co nejmenší. Výstupem má být nějaké optimální rozložení operací pro
jednotlivé stroje (např. ``[[3,3,2],[6,4]]``, což znamená, že použijeme dva
stroje, první z nich vykoná operace trvající ``3``, ``3`` a ``2`` úseky, druhý
operace trvající ``6`` a ``4`` časové úseky, obojí se vejde do limitu ``10``
časových úseků / stroj).

Řešení:

```prolog
sum(List, Sum) :-
    sum_(List, 0, Sum).

sum_([], Acc, Acc).
sum_([X | Xs], Acc, Out) :-
    NewAcc is Acc + X,
    sum_(Xs, NewAcc, Out).

plan(Times, MaxTime, Plan) :-
    msort(Times, TimesSorted),
    plan_(TimesSorted, MaxTime, [], Plan).

plan_([], _, Plan, Plan).
plan_([T | Ts], MaxTime, Plan, Out) :-
    extend_plan(Plan, MaxTime, T, NewPlan),
    plan_(Ts, MaxTime, NewPlan, Out).

extend_plan([], MaxTime, T, [[T]]) :-
    MaxTime >= T.
extend_plan([ P | Ps ], MaxTime, T, Out) :-
    sum(P, PlanTime),
    Free is MaxTime - PlanTime,
    (
        Free >= T
    ->
        Out = [ [ T | P ] | Ps ]
    ;
        extend_plan(Ps, MaxTime, T, Ans),
        Out = [ P | Ans ]
    ).
```

### Listy stromu podle počtu kroků vpravo

Zdroj: [MFF Forum: Zkouška 13. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10963)

Máte zadaný binární strom (klasická ``tree(vlevo, hodnota, vpravo)`` notace).
Roztřiďte vrcholy podle toho, kolikrát musíme jít doprava, než je objevíme.

Řešení:

```prolog
group_by_right_turns(null, []).
group_by_right_turns(tree(Left, Val, Right), Out) :-
    group_by_right_turns(Left, LeftAns),
    group_by_right_turns(Right, RightAns),

    merge_(LeftAns, [ [Val] | RightAns], Out).

merge_(Xs, [], Xs).
merge_([], Ys, Ys).
merge_([X | Xs], [ Y | Ys], [ Z | Ans ]) :-
    append(X, Y, Z),
    merge_(Xs, Ys, Ans).
```

### Maximální párování dle inkluze

Zdroj: [MFF Forum: Zkouška 28.6.2016 - Dvořák, Hric](http://forum.matfyz.info/viewtopic.php?f=169&t=10993)

Napište predikát ``parovani(+G, +H, -P)``, který bere neorientovaný graf ``G``
bez smyček (tj. reflexivních hran) zadaný jako seznam následníků, hranu ``H`` v
podobě ``(v1-v2)`` a vydá co do inkluze maximální párování obsahující zadanou
hranu ``H`` (pozor: nikoli největší párování, ale pouze maximální co do
inkluze).

Například:

```prolog
?- parovani([a-[b,c,d],b-[a,c],c-[a,b,d],d-[a,c],e-[]],a-d,P)
P = [a-d,b-c].
```

Řešení:

```prolog
parovani(Graph, Edge, MaxMatching) :-
    Edge = U-V,
    select(U-_, Graph, Tmp),
    select(V-_, Tmp, RestOfGraph),

    parovani_(RestOfGraph, [U, V], Ans),
    MaxMatching = [ Edge | Ans ].

parovani_([], _, []).
parovani_(Graph, Taken, [ U-V | Ans ]) :-
    select(U-Neighbours, Graph, RestOfGraph),
    member(V, Neighbours),
    \+ member(V, Taken),
    parovani_(RestOfGraph, [U, V | Taken], Ans),
    !.
parovani_(_, _, []).
```

### Generování všech možných výrazů

Zdroj: [MFF Forum: Zkouška 30. 05. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10934)

Na vstupu dostaneme posloupnost čísel a číslo ``V``. Máme vrátit všechny možné
matematické výrazy, které lze z dané posloupnosti postavit pomocí operátorů
``+``, ``-``, ``*``, ``//`` a ``závorek``, a jejichž hodnota je ``V``. Výraz
musí využít všechna zadaná čísla, a jejich pořadí nesmí měnit. Dále si máme
dávat pozor, abychom ve výrazu nedělili nulou.

Řešení:

```prolog
gen_expr(List, V, Expr) :-
    gen_expr_(List, Expr),
    V is Expr.

gen_expr_([Expr], Expr).
gen_expr_(Xs, Expr) :-
    select(X, Xs, Ys),
    select(Y, Ys, Zs),
    !,
    (
        E = X + Y
    ;
        E = X - Y
    ;
        E = X * Y
    ;
        Denom is Y, Denom \= 0, E = X // Y
    ),
    gen_expr_([E | Zs], Expr).
```

### !! Zlepšení řezu

Zdroj: [MFF Forum: Zkouška 19.06.2015 - Dvořák, Hric](http://forum.matfyz.info/viewtopic.php?f=169&t=10536)

Napište predikát ``zlepsirez(+Graf, +Vrcholy1, +Vrcholy2, -OutV)``, který pro
zadaný ohodnocený neorientovaný graf ``Graf`` a řez (definovaný pomocí dvou
disjunktních množin vrcholů ``Vrcholy1`` a ``Vrcholy2``) najde vrchol, který
když přesuneme do opačné skupiny vrcholů řezu, tak dostaneme řez s lepší cenou.

### Ohodnocení stromu post- a pre-order

Zdroj: [MFF Forum: Zkouška 2. 6. 2015 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10479)

Definujte predikát ``transverse(+Strom,-OhodnocenýStrom)``, který zkopíruje
strukturu stromu ``Strom`` do ``OhodnocenýStrom`` s tím, že ke každému vrcholu
přidá číslo ``N``, které znamená pořadí v preOrder průchodu a číslo ``M``,
které znamená pořadí v postOrder průchodu. Ideálně jedním průchodem stromem.

Příklad

```prolog
?- transverse(t(t(nil,l,nil),v,t(nil,p,nil)),X).
X = t(t(nil,l-2-1,nil),v-1-3,t(nil,p-3-2,nil))
```

Řešení:

```prolog
transverse(Tree, Out) :-
    transverse(Tree, 0, 0, _, _, Out),
    !.

transverse(nil, PreOrder, PostOrder, PreOrder, PostOrder, nil).
transverse(Tree, PreOrderIn, PostOrderIn, PreOrderOut, PostOrderOut, Out) :-
    Tree = t(Left, Val, Right),
    NewPreOrder is PreOrderIn + 1,
    transverse(
        Left,
        NewPreOrder,
        PostOrderIn,
        PreOrderOutLeft,
        PostOrderOutLeft,
        LeftAns
    ),
    transverse(
        Right,
        PreOrderOutLeft,
        PostOrderOutLeft,
        PreOrderOutRight,
        PostOrderOutRight,
        RightAns
    ),

    PreOrderOut = PreOrderOutRight,
    PostOrderOut is PostOrderOutRight + 1,
    Out = t(LeftAns, Val-NewPreOrder-PostOrderOut, RightAns).
```

### !! Rotace seznamu

Zdroj: [MFF Forum: Zkouška 25. 5. 2014 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10462)

- napište predikát ``rotace/2`` na rotování seznamu. Nesmíte použít žádné
  pomocné predikáty. (V lineárním čase) (pouze 3 verze)
- v konstantním čase, jakou potřebujete strukturu? Ukažte na ``[1,2,3]``
- napište ``rotace/2`` (pouze 2 verze) v konstantním čase

### Izomorfizmus bin. stromů s popisem

Zdroj: [MFF Forum: Zkouška 12.7.2021](http://forum.matfyz.info/viewtopic.php?f=169&t=12193)

Jsou zadány dva binární (zakořeněné) stromy ``S`` a ``T`` s ohodnocenými
vrcholy, přičemž ohodnocení vrcholů se může opakovat. Definujte predikát
``iso/3``, který zjistí, zdali jsou tyto stromy isomorfní a vydá popis
transformace. Volání je ``iso(+S,+T, -Popis)``, kde ve třetím argumentu bude
popis. Popis je strom stejného tvaru jako ``S`` a ve vrcholech má boolovské
hodnoty ``true`` a ``false``. Hodnota ``true`` ve vrcholu znamená, že se děti
vrcholu v ``S`` mají přehodit, abychom dostali ``T``.

Dva binární stromy jsou isomofní, pokud lze jeden získat z druhého permutací
dětí libovolných vrcholů stromu, tj. vyměněním nebo nevyměněním podstromů
vrcholu.

1. Navrhněte reprezentaci binárního (zakořeněného) stromu s ohodnocenými
   vrcholy v jazyce Prolog. Vaši reprezentaci ukažte na příkladě.
2. Definujte predikát ``iso/3``.
3. Je některý z predikátů, které ve vašem řešení používáte (ať už vámi
   definovaných či knihovních), nedeterministický? Je predikát ``iso/3``
   nedeterministický? Lze ho zdeterminičtit (a jak?), pokud nám stačí nejvýš
   jedno řešení?

Příklad:

```none
  S= d                 T= d                Popis= t
   /---\                /---\                   /---\
  b     e              e     b                 f     t
 / \   / \            / \   / \               / \   / \
a   c f   g          g   f a   c             f   f f   f
```

```prolog
S = t(
    t(
        t(nil, a, nil),
        b,
        t(nil, c, nil)
    ),
    d,
    t(
        t(nil, f, nil),
        e,
        t(nil, g, nil)
    )
).

T = t(
    t(
        t(nil, g, nil),
        e,
        t(nil, f, nil)
    ),
    d,
    t(
        t(nil, a, nil),
        b,
        t(nil, c, nil)
    )
).
```

Řešení:

```prolog
iso(TreeA, TreeB, Transform) :-
    transform(TreeA, Transform, TreeB),
    !.

transform(nil, nil, nil).
transform(t(Left, Val, Right), t(TransLeft, Bool, TransRight), TreeOut) :-
    transform(Left, TransLeft, LeftOut),
    transform(Right, TransRight, RightOut),
    (
        Bool = false, TreeOut = t(LeftOut, Val, RightOut)
    ;
        Bool = true, TreeOut = t(RightOut, Val, LeftOut)
    ).
```

### FirstFit

Dostanete informaci o obsazené paměti jako seznam dvojic ``zacatek-konec`` o
jednotlivých obsazených úsecích. Úseky jsou v seznamu uspořádány vzestupné a
nenavazují bezprostředně na sebe (tj. navazující úseky se spojí) a tyto
invarianty chcete udržovat.

Dále dostanete seznam délek úseků, které máte naalokovat.

Napište predikát

```prolog
firstFit(+Aalokovat, +Obsazeno, -Umisteni, -ObsszenoO)
```

, který naalokuje postupně všechny požadavky z ``Alokovat`` metodou firstFit,
tj. alokuje na první místo, kde se úsek vejde a tím ho obsadí. Vydejte nový
popis obsazených úseků, ve tvaru splňujicím invariant, a popis umístění jako
seznam dvojic ``delkaUseku-umisteni`` ve stejném pořadíjako v ``Alokova``.

Příklad:

```prolog
?- firstFit([100,117,501, 10-50, 1P0-150, 250-1001, U, O).
U = [100-150,10-50,50-100],
O = [0-60, 100-150]
```

Řešení:

```prolog
first_fit([], Obsazeno, [], Obsazeno).
first_fit([H | T], Obsazeno, [ H-U | UmistnenoAns], ObsazenoOut) :-
    first_fit_one(H, Obsazeno, 0-0, U, ObsazenoTmp),
    first_fit(T, ObsazenoTmp, UmistnenoAns, ObsazenoOut),
    !.

first_fit_one(Size, [], LastFrom-LastTo, LastTo, [LastFrom-NewTo]) :-
    NewTo is LastTo + Size.
first_fit_one(Size, [From-To | Rest], LastFrom-LastTo, LastTo, ObsazenoOut) :-
    Free is From - LastTo,
    Free >= Size,

    ObsazenoOut = ObsazenoOut_,
    (
        Free is Size
    ->
        ObsazenoOut_ = [ LastFrom-To | Rest]
    ;
        NewTo is LastTo + Size,
        ObsazenoOut_ = [ LastFrom-NewTo, From-To | Rest]
    ).
first_fit_one(Size, [From-To | Rest], LastFrom-LastTo, OutPos, ObsazenoOut) :-
    Free is From - LastTo,
    Free < Size,

    first_fit_one(Size, Rest, From-To, OutPos, ObsazenoAns),

    ObsazenoOut = ObsazenoOut_,
    (
        0 is LastTo
    ->
        ObsazenoOut_ = ObsazenoAns
    ;
        ObsazenoOut_ = [LastFrom-LastTo | ObsazenoAns]
    ).
```

### Otočení v sekvenci

Na vstupu je daný seznam ``S`` nějakých položek, například RNA bází. Chcete
vydat seznam seznamú položek ``Vs`` jako seznam výsledků, který vznikne
otočením nějaké souvislé části ``S`` délky aspoň ``2`` všemi možnými zpüsoby.
Napište predikát ``otoceni(+S, -Vs)``.

Přiklad:

```Prolog
?- otoceni([ a, c, g, t], Vs).
Vs = [[c, a, g, t], [g, c, a, t], [t, g, c, a], [a, t, g, c], [a, c, t, g]]
```

Řešení:

```prolog
je_otoceni(List, Out) :-
    append(Front, MidBack, List),
    append(Mid, Back, MidBack),

    length(Mid, N),
    N >= 2,

    reverse(Mid, MidRev),

    append(MidRev, Back, Tmp),
    append(Front, Tmp, Out).


otoceni(List, Out) :-
    otoceni_(List, [], Out).

otoceni_(Xs, Acc, Out) :-
    je_otoceni(Xs, P),
    \+ member(P, Acc),
    otoceni_(Xs, [P | Acc], Out),
    !.
otoceni_(_, Acc, Acc) :- !.
```

## Haskell

### Největší součet souvislé podposloupnosti

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095)

1. Pro zadanou posloupnost čísel najděte spojitý úsek, jehož součet je
   největší. Vydejte souřadnice začátku a konce úseku a dosažený součet.

    ```haskell
    soucty :: Num a => [a] → (Int, Int, a)
    ```

    Pokuste se o nějakou optimalizaci, tj. nepočítejte součty hrubou silou (zcela
    samostatně).

    Příklad: (indexováno od 0)

    ```haskell
    > soucty [-1,1,2,3,-4]
     (1,3,6)
    ```

2. Jaký význam má část ``Num a =>`` v definici funkce soucty ? Proč tam musí
   být?
3. Uveďte dvě možné konkrétní hodnoty proměnné a z typu funkce soucty.
4. Lze definovat ``Num a``taky pro uživatelské typy nebo musíme použít pouze
   předdefinované/vestavěné? Lze naši funkci soucty použít pro nějaký
   uživatelský typ na místě ``a`` ? (Proč ano/ne?)

Řešení:

```haskell
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

Zdroj: [MFF Forum: Zkoužka 4.9.2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12095)

Binární vyhledávací strom je zadán jako seznam hodnot vrcholů v pořadí
preorder. Definujte funkci

```haskell
readBt :: Ord a => [a] -> Bt a
```

která ze zadaného seznamu zrekonstruuje původní strom typu

```haskell
data Bt a = Void
          | Node (Bt a) a (Bt a)
```

Připomeňme, že v binárním vyhledávacím stromu platí pro každý vrchol ``v``, že
všechny hodnoty v levém, resp. pravém podstromu v jsou menší, resp. větší nežli
hodnota ``v``. Odtud plyne, že původní strom je zadaným seznamem určen
jednoznačně.

Příklad:

```haskell
> readBt [5, 2, 4, 9]
Node (Node Void 2 (Node Void 4 Void)) 5 (Node Void 9 Void)
```

1. Definujte funkci ``readBt``.
2. Je ve vašem řešení použita nějaká funkce vyššího řádu (funkce s
   funkcionálními argumenty)? Pokud ne, dala by se zde nějaká smysluplně
   použít?
3. Je ve vašem řešení použita notace stručných seznamů (list comprehension),
   tj. ``[... | ...]`` ? Pokud ne, dala by se zde smysluplně použít?

Řešení:

```haskell
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

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089)

Cílem je definovat funkci ``ascending,`` která na vstupu obdrží seznam hodnot
(libovolného typu) a vrátí zpět seznam posloupností, který splňuje:

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
3. Jak byste zobecnili tuto funkci tak, aby ji bylo možné použít s libovolným
   porovnávacím operátorem?
4. Bude vaše definice fungovat i na nekonečných seznamech? Pokud ano,
   vysvětlete proč. Pokud ne, dala by se vaše definice takto upravit?
   Zdůvodněte proč.

Řešení:

```haskell
ascending :: Ord a => [a] -> [[a]]
ascending [] = []
ascending xs = ys : ascending zs where (ys, zs) = takeAscending (<) xs

takeAscending :: Ord a => (a -> a -> Bool) -> [a] -> ([a], [a])
takeAscending cmp []  = ([], [])
takeAscending cmp [x] = ([x], [])
takeAscending cmp (x1 : x2 : xs) | cmp x1 x2   = (x1 : asc, rest)
                               | otherwise = ([x1], x2 : xs)
  where (asc, rest) = takeAscending cmp (x2 : xs)
```

### Stromové operace

Zdroj: [MFF Forum: Zkouška 16. 7. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12089)

 1. Definujte datový typ pro binární stromy.
    - Hodnoty jsou uloženy ve vnitřních uzlech.
    - Pokuste se o co nejobecnější definici.
    - Nezapomeňte na reprezentaci prázdného stromu.
 2. Definujte funkci ``replicateT``. Výsledkem ``replicateT n a`` je binární
    strom, který obsahuje ``n`` kopií hodnoty ``a``.
    - Výsledný strom by měl mít minimální možnou hloubku. Např. strom
      ``replicateT 7 a`` by měl mít hloubku 3.
 3. Definujte funkci ``zipWithT`` jako zobecnění funkce ``zipWith``. ``zipWithT
    f t1 t2`` sloučí prvky stromů ``t1`` a ``t2`` na stejných pozicích pomocí
    funkce f.
    - Pokud nemá nějaký prvek z jednoho stromu odpovídající prvek na stejné
      pozici v druhém stromě, tak jej do výsledného stromu nepřidávejte. Např.
      pro prázdný strom empty by mělo platit ``zipWithT f t empty == empty`` a
      ``zipWithT f empty t == empty``.
 4. Pomocí ``replicateT`` a ``zipWithT`` definujte funkci ``cut``. Funkce ``cut
    n t`` odstraní ze stromu ``t`` všechny vrcholy, jejichž hloubka je ostře
    větší než ``n``.

Řešení:

```haskell
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

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078)

Cílem je definovat binární funkci klouzave, která

- obdrží na vstupu posloupnost čísel a přirozené číslo ``n``
- a vrátí posloupnost klouzavých průměrů řádu ``n``, tj. aritmetických průměrů
  ``n`` sousedních prvků.

Příklad:

```haskell
klouzave [1.5, 2.5, 3.5, 4.5, 5.5] 3
[2.5,3.5,4.5]
```

1. Definujte typovou signaturu funkce ``klouzave``
2. Definujte vlastní funkci s explicitním využitím rekurze
3. Sestavte alternativní definici, tentokráte bez explicitního použití rekurze,
   přitom můžete využívat libovolné knihovní funkce z přiloženého seznamu.
4. Vyhýbá se alespoň jedna z vašich definic opakovaným výpočtům? Pokud ne, dala
   by se takto upravit? Zdůvodněte.
5. Bude některá z vašich definic fungovat i na nekonečných seznamech? Pokud
   ano, vysvětlete proč. Pokud ne, dala by se některá z vašich definic takto
   upravit? Zdůvodněte.

```haskell
take 5 $ klouzave [1..] 10
[5.5,6.5,7.5,8.5,9.5]
```

Řešení:

```haskell
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

Zdroj: [MFF Forum: Zkouška 24. 6. 2020](http://forum.matfyz.info/viewtopic.php?f=169&t=12078)

Cílem toho problému je zobecnit funkce ``foldr`` / ``foldl`` na obecné kořenové stromy.

1. Definujte datový typ pro reprezentaci obecných kořenových stromů s
   ohodnocenými vrcholy:
    - snažte se o co nejobecnější definici
    - nezapomeňte na reprezentaci prázdného stromu
2. Funkce ``foldl`` a ``foldr`` zobecněte na funkci ``foldT``, která bude -
   namísto seznamu - procházet stromem ve vaší reprezentaci popsané v **1.**.
3. Pomocí funkce fold definujte funkci ``arita,`` která vrátí ``aritu`` (tj.
   maximální počet dětí přes všechny vrcholy) zadaného kořenového stromu.
4. Pomocí funkce ``foldT`` definujte funkci ``pdc``, která vrátí průměrnou
   délku cesty z kořene do listu (tj. součet délek všech cest z kořene do listu
   / počet listů).

Řešení:

```haskell
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

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066)

Rozdělte zadaný binární vyhledávací strom ``T`` na ``n+1`` binárních
vyhledávacích stromů ``T_0, .. , T_n`` podle zadaných vstupních hodnot ``k_i``,
``1 <= i <= n`` tak, že ve stromě ``T_i`` jsou hodnoty ``x``, ``k_i <= x <
k_i+1``, pokud jsou nerovnosti aplikovatelné.

Obrázek
![priklad vstupu](https://i.imgur.com/sNvonyI.png)

Snažte se o efektivitu, celé podstromy patřící do jednoho pruhu zpracujte najednou.

1. Definujte datový typ pro reprezentaci binárních vyhledávacích stromů. Snažte
   se o co nejobecnější definici.
2. Definujte typovou signaturu funkce ``pruhy``, včetně typových tříd.
3. Funkci ``pruhy`` definujte. Budete-li používat pomocné funkce, u každé
   popište její význam.
4. Pokuste se stručně zdůvodnit korektnost vaší defnice.

Řešení:

```haskell
data BTree a = Nil
          | BTree (BTree a) a (BTree a)
          deriving (Eq, Show)

cutUpTo :: Ord a => a -> BTree a -> (BTree a, BTree a)
cutUpTo _ Nil = (Nil, Nil)
cutUpTo max (BTree left val right)
  | val >= max
  = let (ans, restAns) = cutUpTo max left in (ans, BTree restAns val right)
  | otherwise
  = let (ans, restAns) = cutUpTo max right in (BTree left val ans, restAns)

pruhy :: Ord a => [a] -> BTree a -> [BTree a]
pruhy [] tree = [tree]
pruhy (x : xs) tree =
  let (part, rest) = cutUpTo x tree in part : pruhy xs rest


testTree =
  BTree (BTree (BTree Nil 1 Nil) 2 (BTree Nil 4 Nil)) 5 (BTree Nil 6 Nil)
```

### Run-length encoding/decoding

Zdroj: [MFF Forum: Zkouška 12. 6. 2020 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=12066)

Definujte funkce ``rle`` a ``rld``, které realizují run-length encoding a
decoding. Funkce

```haskell
rle :: Eq a => [a] -> [Either a (a,Int)]
```

zakóduje co nejdelší úseky stejných prvků ve vstupním seznamu do dvojice
``(prvek, počet)`` typu ``Either`` s datovým konstruktorem ``Right``.

Pokud je prvek v úseku sám, kóduje se pouze prvek vnořený do typu ``Either`` s
datovým konstruktorem ``Left``.

Příklad:

```haskell
> rle ”abbcccda”
[Left 'a', Right ('b',2), Right ('c',3), Left 'd', Left 'a']
```

1. Definujte funkci ``rle`` s využitím rekurze, ale bez použití stručných
   seznamů či funkcí vyšších řádů (funkce s funkcionálními parametry).
2. Definujte funkci ``rle`` bez explicitního využití rekurze, ale za použití
   stručných seznamů či funkcí vyšších řádů.
3. Definujte typovou signaturu funkce ``rld``, která realizuje dekompresi, tj.
   převod ze seznamu úseků na původní seznam prvků.
4. Definujte funkci ``rld.`` Použijte přitom funkci ``map`` či ``concat``.
5. Bude některá z funkcí fungovat i na nekonečných seznamech? Proč ano nebo
   proč ne?

Řešení:

```haskell
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

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977)

Definujte funkce:

```none
prevod1 cislo puvodni
```

pro převod čísla z číselné soustavy o základu ``puvodni`` do ``dekadické
číselné soustavy``, a

```none
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

    ```haskell
    prevod1 ::
    prevod2 ::
    ```

2. Definujte funkci ``prevod1`` s využitím rekurze.
3. Sestavte alternativní definici funkce prevod1 s využitím alespoň jedné z
   funkcí ``map``, ``filter``, ``foldr`` či ``foldl``, ale bez (explicitního)
   použití rekurze.
4. Definujte funkci ``prevod2`` s využitím funkce ``unfold`` definované
   následovně:

    ```haskell
    unfold :: (t -> Bool) -> (t -> (a, t)) -> t -> [a]
    unfold done step x =  if done x then []
                                    else let (y,ys) = step x
                                         in y: unfold done step ys
    ```

Řešení:

```haskell
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

Zdroj: [MFF Forum: 10. 9. 2019 - Dvořák](http://forum.matfyz.info/viewtopic.php?f=169&t=11977)

Definujte unární funkci rady, která obdrží multiplikativní tabulku grupy jako
matici prvků. První řádek matice obsahuje násobení grupovou jednotkou ``e`` a
pořadí prvků odpovídající řádkům a sloupcům je stejné. Vydá seznam všech prvků
spolu s jejich řády.

Řád prvku ``p`` je nejmenší přirozené číslo ``n`` takové, že ``n``-tá mocnina
``p`` je rovna ``e``.

1. Definujte typovou signaturu funkce ``rady``.
2. Funkci ``rady`` definujte.

Příklad:

```haskell
> rady [["e","a","b"], ["a","b","e"], ["b","e","a"]]
[("e",1), ("a",3), ("b",3)]
```

Řešení:

```haskell
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

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969)

Je dána číselná matice ``A``. Definujte funkci

```haskell
kumulace :: Num a => [[a]] -> [[a]]
```

která z matice ``A`` vyrobí matici ``B`` stejných rozměrů (viz příklad níže).

Každý prvek na souřadnicích ``(i,j)`` bude roven součtu všech hodnot v
submatici s levým horním rohem ``(0,0)``a pravým dolním rohem ``(i,j)``.

*Poznámka: Snažte se vyhnout opakování stejných výpočtů.*

Příklad:

```haskell
> kumulace[[1,1,1],[1,2,1],[0,1,0],[1,1,-4]]
[[1,2,3],[2,5,7],[2,6,8],[3,8,6]]
```

Řešení:

```haskell
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

Zdroj: [MFF Forum: 24.6.2019 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11969)

Hypergraf je zadán množinou vrcholů a množinou hyperhran, což jsou alespoň
dvouprvkové podmnožiny množiny vrcholů. Naší cílem je definovat funkci
doplnění, která doplní do hypergrafu ``H`` všechny dvouprvkové (hyper)hrany pro
ty dvojice vrcholů, které nejsou společně obsaženy v žádné hyperhraně vstupního
hypergrafu ``H``. Funkce tedy např. z hypergrafu s vrcholy ``{1,2,3,4,5}`` a
hyperhranani ``{1,3,5}`` a ``{2,3,4}`` vytvoří hypergraf se stejnými vrcholy a
hyperhranami ``{1,3,5},{2,3,4},{1,2},{1,4},{5,2}`` a ``{5,4}``

1. Definujte datový typ pro reprezentaci hypergrafu. Pokuste se o co
   nejobecnější definici (vrcholy mohou být reprezentovány nejen čísly, ale i
   znaky, řetězci apod.)
2. Specifikujte typovou signaturu funkce

    ```haskell
    doplneni ::
    ```

3. Funkci definujte.

Řešení:

```haskell
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

### Analýza textu

Zdroj: [Zkouška 10.6.2019 (Dvořák + Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=11954)

Na vstupu je zadán text jako hodnota typu ``String``. Naším cílem je definovat
binární funkci ``stat text n``, která

- obdrží takový text a přirozené číslo ``n``
- vrátí všechna slova z tohoto textu o délce alespoň ``n``, setříděná lexikograficky
- každé slovo s čísly řádků, kde se slovo vyskytuje

Řádky jsou ukončeny znakem ``'\n'``. Slovo je každý maximální podřetězec textu
neobsahující mezeru ``' '``, tabulátor ``'\t'`` či konec řádku ``'\n'``.

1. Definujte datovou strukturu pro reprezentaci oboru hodnot funkce stat
   (pomocí data nebo type).
2. Definujte typovou signaturu funkce stat s použití datové struktury z **1.**.
3. Funkci stat definujte.

Řešení:

```haskell
newtype Stat = Stat [(Int, String)]
  deriving (Eq, Show)

lines' :: String -> [String]
lines' "" = []
lines' ss = line : lines' rest
 where
  line = takeWhile (/= '\n') ss
  rest = dropWhile (/= '\n') ss

words' :: String -> [String]
words' "" = []
words' ss = case takeWhile (not . isSpace) ss of
              [] -> []
              word -> word : words' rest
 where
  isSpace = flip elem [' ', '\t', '\n']

  tmp     = dropWhile (not . flip elem [' ', '\t', '\n']) ss
  rest    = dropWhile isSpace tmp

sortBy :: Ord b => (a -> b) -> [a] -> [a]
sortBy _ []           = []
sortBy f (pivot : xs) = left ++ [pivot] ++ right
 where
  left  = filter (\y -> f y <= f pivot) xs
  right = filter (\y -> f y > f pivot) xs

stat :: String -> Int -> Stat
stat text n = Stat sortedWords
 where
  numberedLines = zip [1 ..] $ lines' text
  linesToWords  = concatMap
    (\(line_no, line) -> zip (repeat line_no) (words' line))
    numberedLines
  filteredWords = filter (\(_, word) -> length word >= n) linesToWords
  sortedWords   = sortBy snd filteredWords
```

### Označkování stromu

Zdroj: [MFF Forum: Zkouška 21.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11747)

Ohodnotit vrcholy obecného n-arní stromu v pořadí v jakém jsou vrcholy
uzavírány, takže post-fixově. Měla se napsat datová struktura pro strom, typová
hlavička fce a pak tu funkci implementovat:

```haskell
data Tree a = Nil | Tree a [Tree a]
label :: Tree a -> Tree (a, Int)
```

Příklad

```haskell
> label $ Tree 1 [Tree 1 [Nil],Nil,Tree 2 [Tree 4 [Tree 5 [Nil],Tree 6 [Nil]]],Tree 3 [Nil,Nil,Nil,Nil]]
Tree (1,7) [Tree (1,1) [Nil],Nil,Tree (2,5) [Tree (4,4) [Tree (5,2) [Nil],Tree (6,3) [Nil]]],Tree (3,6) [Nil,Nil,Nil,Nil]]
```

Řešení:

```haskell
data Tree a = Nil | Tree a [Tree a]
  deriving (Eq, Show)

label :: Tree a -> Tree (a, Int)
label = snd . label' 0

label' :: Int -> Tree a -> (Int, Tree (a, Int))
label' n Nil           = (n, Nil)
label' n (Tree val ts) = (newN, Tree (val, newN) ansTs)
 where
  (ansN, ansTs) = sequentialLabel n ts
  newN          = ansN + 1

sequentialLabel :: Int -> [Tree a] -> (Int, [Tree (a, Int)])
sequentialLabel n [] = (n, [])
sequentialLabel n (Nil : ts) =
  let (ansN, ansTs) = sequentialLabel n ts in (ansN, Nil : ansTs)
sequentialLabel n (tree : ts) =
  let (ansN1, ansT ) = label' n tree
      (ansN2, ansTs) = sequentialLabel ansN1 ts
  in  (ansN2, ansT : ansTs)
```

### Změna některých prvků

Zdroj: [MFF Forum: Zkouška 26.6.2018](http://forum.matfyz.info/viewtopic.php?f=169&t=11756)

Definujte funkci ``change``, která:

- obdrží seznam, který reprezentuje abecedu povolených prvků, které lze využít
  k modifikaci
- obdrží seznam ``xs`` pro modifikaci
- vrátí seznam všech modifikací vstupního seznamu ``xs``, které se od něho liší
  v právě 3 prvcích.

Příklad:

```haskell
> change3 "ab" "aabe"
["bbae", "bbba", "bbbb", "baaa", "baab", "abaa", "abab"]
```

1. Definujte typ funkce ``change3`` co nejobecněji (včetně případných typových tříd)
2. Definujte funkci ``change3``

Řešení:

```haskell
change3 :: Eq a => [a] -> [a] -> [[a]]
change3 cs xs = map snd $ filter (\(count, _) -> count == 3) $ change' cs xs

change' :: Eq a => [a] -> [a] -> [(Int, [a])]
change' _  []       = [(0, [])]
change' cs (x : xs) = changed ++ notChanged
 where
  ans        = change' cs xs

  other      = filter (/= x) cs
  notChanged = map (\(count, ys) -> (count, x : ys)) ans
  changed =
    map (\(c, (countChanged, ys)) -> (countChanged + 1, c : ys))
      $ [ (c, p) | c <- other, p <- ans ]
```

### Největší kladná podmatice

Zdroj: [MFF Forum: Zkouška 6. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10961)

Máme zadanou matici (jako seznam seznamů). Naším cílem je vypsat seznam všech
dvojic ``(x, y)`` takových, že podmatice ``(1, 1):(x, y)`` bude obsahovat pouze
kladné hodnoty. Dvojice ``(x, y)`` musí být vždy nejvyšší možné (t. j. nelze
ani v jedne souradnici zvětšit)

Řešení:

```haskell
scan :: (b -> a -> b) -> b -> [a] -> [b]
scan _ acc []       = [acc]
scan f acc (x : xs) = acc : scan f (f acc x) xs

maxPositive :: (Num a, Ord a) => [[a]] -> [(Int, Int)]
maxPositive []   = []
maxPositive [[]] = []
maxPositive matrix =
  [ (row, count)
  | (row, count, jump) <- zip3 [1 ..] maxCounts jumpDown
  , count > 0
  , jump
  ]
 where
  posCounts = map (length . takeWhile (> 0)) matrix
  maxCounts = drop 1 $ scan min (maxBound :: Int) posCounts
  jumpDown =
    [True]
      ++ [ curr > next
         | (curr, next) <- zip (drop 1 posCounts) (drop 2 posCounts)
         ]
      ++ [True]
```

### Stromový fold

Zdroj: [MFF Forum: Zkouška 6. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10961)

1. Napiste fold pro binarni stromy

    ```haskell
    data Tree a = Leaf a | Tree (Tree a) (Tree a)
    fold :: (a -> b) -> (b -> b -> b) -> Tree a -> b
    ```

2. Napiste one-liner funkci, ktera vypise minimum a maximum z celeho stromu
   pomoci vami napsaneho foldu.
3. Napiste hlavicku funkce z **2.**

Řešení:

```haskell
data Tree a = Leaf a | Tree (Tree a) (Tree a)
  deriving (Eq, Show)

fold :: (a -> b) -> (b -> b -> b) -> Tree a -> b
fold f _    (Leaf a         ) = f a
fold f comb (Tree left right) = comb leftAns rightAns
 where
  leftAns  = fold f comb left
  rightAns = fold f comb right

minmaxT :: Ord a => Tree a -> (a, a)
minmaxT = fold
  (\a -> (a, a))
  (\(minLeft, maxLeft) (minRight, maxRight) ->
    (min minLeft minRight, max maxLeft maxRight)
  )
```

### Tetris

Zdroj: [MFF Forum: Zkouška 22.6.](http://forum.matfyz.info/viewtopic.php?f=169&t=11412)

Máme obdélníkovou tabulku uloženou po řádcích jako seznam seznamů Intů. Vymažte
z ní všechny sloupce, které neobsahují žádnou nulu.

Řešení:

```haskell
transpose :: [[a]] -> [[a]]
transpose []     = []
transpose matrix = case concatMap (take 1) matrix of
  []  -> []
  col -> col : transpose (map (drop 1) matrix)

tetris :: (Eq a, Num a) => [[a]] -> [[a]]
tetris = transpose . removeFull . transpose
  where removeFull rows = [ row | row <- rows, 0 `elem` row ]
```

### Splnění všech podmínek

Zdroj: [MFF Forum: Zkouska 20.9.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11466)

Funkce ``podminky`` dostává seznam podmínek jedné proměnné a seznam hodnot.
Vydává seznam seznamů hodnot, kde ``i``-tý seznam na výstupu obsahuje hodnoty,
pro které byla splněna ``i``-tá podmínka a nebyly splněny předchozí podmínky.
Hodnoty, pro které nebyla splněna žádná podmínka, se zahodí.

Příklad:

```haskell
> podminky [even,(>5),(==3)] [0..9]
[[0,2,4,6,8],[7,9],[3]]
```

1. Napište typovou signaturu funkce podmínky (co nejobecnější, včetně
   případných typových tříd).
2. Napište funkci ``podminky``.

Řešení:

```haskell
podminky :: [a -> Bool] -> [a] -> [[a]]
podminky []       _  = []
podminky (f : fs) xs = filter f xs : podminky fs rest
  where rest = filter (not . f) xs
```

### Stromový take

Zdroj: [MFF Forum: Zkouska 20.9.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11466)

Cílem tohoto problému je zobecnit standardní funkci ``take`` na funkci
``takeTree``, která

- obdrží obecný kořenový strom a dvě přirozená čísla ``n`` a ``m``
- odstraní ve stromě všechny vrcholy ve hloubce větší než ``m`` (hloubka
  vrcholu ``v`` je počet hran na cestě z kořene do ``v``)
- pro každý vrchol, který má více než ``n`` dětí, odstraní všechny děti (s
  příslušnými podstromy) kromě ``n`` nejlevějších
- výsledný (nejvýše ``n``-ární) strom (hloubky nejvýše ``m``) vrátí.

1. Definujte datový typ pro obecný kořenový strom, v jehož vrcholech jsou
   uloženy prvky typu ``a``.
2. Využijte váš datový typ k definici nekonečného stromu, tj. takového stromu,
   že pro každé přirozené číslo ``i`` buďto existuje vrchol s alespoň ``i``
   dětmi, nebo existuje vrchol ve hloubce alespoň ``i``.
3. Definujte typovou signaturu funkce ``takeTree``.
4. Funkci ``takeTree`` definujte.

Řešení:

```haskell
data Tree a = Nil | Tree a [Tree a]
  deriving (Eq, Show)

infiniteTree :: Tree Int
infiniteTree = go 0 where go n = Tree n (take (n + 1) $ repeat (go (n + 1)))

takeTree :: Int -> Int -> Tree a -> Tree a
takeTree n m = go 0
 where
  go _ Nil = Nil
  go h (Tree val ts)
    | h == m    = Tree val []
    | otherwise = let tsAns = take n $ map (go (h + 1)) ts in Tree val tsAns
```

### Formule

Zdroj: [MFF Forum: Zkouška 13. 9. 2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11457)

Máme typ:

```haskell
data Formule = Konst Bool | Not Formule | And Formule Formule | Or Formule Formule
```

a chceme napsat funci ``gen``, která vygeneruje nekonečný seznam složený z formulí:

Příklad:

```haskell
gen = [ Konst True, Konst Flase, Not True, Not False, And True True, ... ]
```

Řešení:

```haskell
data Formule = Konst Bool
             | Not Formule
             | And Formule Formule
             | Or Formule Formule
             deriving(Eq, Show)

gen :: [Formule]
gen = concat memo

memo = map genN [0 ..]

genN :: Int -> [Formule]
genN 0 = []
genN 1 = [Konst True, Konst False]
genN n = ands ++ ors
 where
  last = memo !! (n - 1)
  nots = [ Not f | f <- last ]
  ands = [ And f s | f <- last, s <- last ]
  ors  = [ Or f s | f <- last, s <- last ]
```

### Převody stromů

Zdroj: [MFF Forum: Zkouška 6. 6. 2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11380)

Máme dva druhy stromů - obecný n-ární:

```haskell
data NTree a = NTree a [NTree a]
```

a n-ární, ve kterém je řečeno, které podstromy jsou vlevo a které vpravo:

```haskell
data UspTree a = UspTree [UspTree a] a [UspTree a]
```

Máme napsat funkci, která obecný n-ární strom převede na uspořádaný strom. V
každém uzlu obecného n-árního stromu na vstupu je kromě hodnoty uložený taky
počet synů, kteří jsou vlevo.

Řešení:

```haskell
data NTree a = NTree a [NTree a]
  deriving (Eq, Show)

data UspTree a = UspTree [UspTree a] a [UspTree a]
  deriving (Eq, Show)

prevodT :: NTree (Int, a) -> UspTree a
prevodT (NTree (n, val) ts) = UspTree (take n ts') val (drop n ts')
  where ts' = map prevodT ts

```

### Počet trojúhelníků

Zdroj: [MFF Forum: Zkouška 29.5.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11357)

- Navrhněte datový typ ``Graf`` a pro reprezentaci konečného neorientovaného
  grafu s vrcholy typu ``a``.
- Definujte funkci ``troj :: Graf a -> Int``, která k takovému grafu vrátí
  počet všech jeho trojúhelníků.

Priklad:

```haskell
> let
testGraph = Graph [0..8]
                  [(0, 1), (0, 3), (1, 0), (1, 2), (1, 3), (2, 1), (2, 4),
                   (3, 0), (3, 1), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4),
                   (5, 6), (6, 5), (6, 7), (6, 8), (7, 6), (7, 8), (8, 6),
                   (8, 7)]
> troj testGraph
[(0,1,3),(6,7,8)]
```

Řešení:

```haskell
type Edge a = (a, a)

data Graph a = Graph [a] [Edge a]
  deriving (Eq, Show)

troj :: Eq a => Graph a -> [(a, a, a)]
troj (Graph xs edges) = troj' edges xs

troj' :: Eq a => [Edge a] -> [a] -> [(a, a, a)]
troj' _ [] = []
troj' edges (a : vertices) =
  [ (a, b, c)
  | n <- [1 .. length vertices - 1]
  , b <- drop (n - 1) $ take n vertices
  , (a, b) `elem` edges
  , c <- drop n vertices
  , (b, c) `elem` edges
  , (c, a) `elem` edges
  ]
  ++ troj' edges vertices
```

### Bag fold

Zdroj: [MFF Forum: Zkouška 29.5.2017](http://forum.matfyz.info/viewtopic.php?f=169&t=11357)

Je dán datový typ

```haskell
data Bag a = Item a | Items [Bag a]
```

1. Definujte funkci ``fold`` pro obecný průchod touto datovou strukturou (to
   ``(a->b)`` tam zastupuje počáteční hodnotu v normálním foldu)

    ```haskell
    fold :: (a -> b) -> ([b] -> b) -> Bag a -> b
    ```

2. Pomocí funkce fold definujte funkci ``listy`` která posbírá všechny hodnoty
   z položek ``Item`` ze všech úrovní zleva.

```haskell
listy :: Bag a -> [a]
```

Příklad:

```haskell
> listy (Items [Item 1,Items [Item 2, Item 3], Items [Items [Item 4]]])
[1,2,3,4]
```

Řešení:

```haskell
data Bag a = Item a | Items [Bag a]
  deriving (Eq, Show)

fold :: (a -> b) -> ([b] -> b) -> Bag a -> b
fold f _    (Item  a ) = f a
fold f comb (Items bs) = comb ans where ans = map (fold f comb) bs

listy :: Bag a -> [a]
listy = fold (: []) concat
```

### Hledání skoku

Zdroj: [MFF Forum: Zkouška 13. 6. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10963)

Máte nějakou funkci, která nabývá jen dvou různých funkčních hodnot. Funkce
přechází někde (nevíme kde) skokově z jedné funkční hodnoty na druhou. Na
vstupu dostanete ``c`` a ``d`` určující ony dvě funkční hodnoty. Dále dostanete
seznam ``(x, y)`` bodů, ve kterých jste funkci změřili s nějakou chybou.

Napište funkci, která na výstupu tyto body rozdělí na levé a pravé (seznam dvou
seznamů) podle toho, které body patří ještě k hodnotě ``c``, a které už k
hodnotě ``d``.

Pozor, je potřeba minimalizovat celkovou odchylku spočtenou jako součet
``(f(x_i) - y_i)^2`` přes všechny body, kde ``f(x)`` je změřená hodnota (ze
seznamu) a ``y`` skutečná hodnota z našeho odhadu.

Řešení:

```haskell
sortBy :: Ord b => (a -> b) -> [a] -> [a]
sortBy _ []       = []
sortBy f (x : xs) = left ++ [x] ++ right
 where
  left  = sortBy f $ filter (\elem -> f elem < f x) xs
  right = sortBy f $ filter (\elem -> f elem >= f x) xs

minimumBy :: Ord b => (a -> b) -> [a] -> a
minimumBy f = (!! 1) . sortBy f

skok :: (Ord a, Num a) => a -> a -> [(a, a)] -> ([(a, a)], [(a, a)])
skok _ _ [] = undefined
skok c d ps = snd $ minimumBy fst cuts
 where
  sorted = sortBy fst ps
  cuts =
    [ (cost left right, (left, right))
    | n <- [0 .. length sorted]
    , let left  = take n sorted
    , let right = drop n sorted
    ]
  cost xs ys =
    sum $ [ (fx - c) ^ 2 | (_, fx) <- xs ] ++ [ (fx - d) ^ 2 | (_, fx) <- ys ]
```

### Násobení/sčítání řídkých polynomů

Zdroj: [MFF Forum: Zkouška 28.6.2016 - Dvořák, Hric](http://forum.matfyz.info/viewtopic.php?f=169&t=10993)

Mějme řídké polynomy reprezentované pomocí ``[(nenulový
koeficient,exponent)]``. Definujte pro ně datový typ (nezapomeňte na nulový
polynom) a napište funkci ``mult`` (i její datovou signaturu), která bude řídké
polynomy násobit.

- *řídký polynom*: u spousty exponentů je nulový koeficient (exponenty prostě
  nejdou po 1, ale skáčou)

```haskell
data Ridky a = Ridky [(Int, a)]
```

Řešení:

```haskell
type Order = Int
type Coeff = Int
newtype Poly = Poly [(Coeff, Order)]
  deriving (Eq, Show)

sortBy :: Ord b => (a -> b) -> [a] -> [a]
sortBy _ []       = []
sortBy f (x : xs) = left ++ [x] ++ right
 where
  left  = sortBy f $ filter (\elem -> f elem < f x) xs
  right = sortBy f $ filter (\elem -> f elem >= f x) xs

groupBy f []           = []
groupBy f xs@(x : xs') = group : groupBy f rest
 where
  group = takeWhile ((== f x) . f) xs
  rest  = dropWhile ((== f x) . f) xs

normalize :: Poly -> Poly
normalize (Poly ps) = Poly p
 where
  groups = groupBy snd $ sortBy snd ps
  p      = [ (sum $ map fst group, snd $ head group) | group <- groups ]

mult :: Poly -> Poly -> Poly
mult (Poly p) (Poly q) = normalize $ Poly ans
  where ans = [ (c1 * c2, o1 + o2) | (c1, o1) <- p, (c2, o2) <- q ]

summ :: Poly -> Poly -> Poly
summ p q = summ' (normalize p) (normalize q)

summ' :: Poly -> Poly -> Poly
summ' (Poly xs) (Poly ys) = Poly $ go xs ys
 where
  go xs [] = xs
  go [] ys = ys
  go xs@(x@(xCoeff, xOrd) : xs') ys@(y@(yCoeff, yOrd) : ys')
    | xOrd == yOrd = (xCoeff + yCoeff, xOrd) : go xs' ys'
    | xOrd > yOrd  = x : go xs' ys
    | otherwise    = y : go xs ys'
```

### Maximo-lexikografické generování všech dvojic

Zdroj: [MFF Forum: Zkouška 30. 05. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10934)

Pro zadané ``K`` máme generovat nekonečný uspořádáný seznam ``K``-tic:
uspořádání je definováno tak, že nejprve se třídí podle maximálního prvku v
daném seznamu. (jakákoliv ``k``-tice, jejíž maximum je menší nebo rovno ``2``
bude před ``k``-ticí obsahující číslo 3). Když mají dvě k-tice stejné maximum,
řadí se lexikograficky.

Příklad pro ``K=2``:

```haskell
[[0,0],[0,1],[1,0],[1,1],[0,2],[1,2],[2,0],[2,1],[2,2],[0,3] ...]
```

Řešení:

```haskell
sort :: Ord a => [a] -> [a]
sort []       = []
sort (x : xs) = left ++ [x] ++ right
 where
  left  = sort $ filter (< x) xs
  right = sort $ filter (>= x) xs

sequencesUpTo :: Int -> Int -> [[Int]]
sequencesUpTo 0   _   = []
sequencesUpTo len max = go len [[]]
 where
  go 0 acc = acc
  go k acc = go (k - 1) [ n : seq | n <- [0 .. max], seq <- acc ]

maxLex :: Int -> [[Int]]
maxLex k = concatMap (sort . withMax) [0 ..]
  where withMax max = [ seq | seq <- sequencesUpTo k max, max `elem` seq ]

```

### Ořezání intervalu z BVS

Zdroj: [MFF Forum: Zkouška 30. 05. 2016 (Dvořák, Hric)](http://forum.matfyz.info/viewtopic.php?f=169&t=10934)

Máme zadaný binární vyhledávací strom (reprezentaci si máme zvolit), a dvě
čísla ``D``, ``H``. Máme vrátit BVS, který vznikl ořezáním vstupního stromu
tak, aby v něm byly pouze hodnoty ``X`` takové, že ``D<=X<=H``.

Řešení:

```haskell
data Tree a = Nil | Tree (Tree a) a (Tree a)
  deriving (Eq, Show)

cut :: (Ord a, Num a) => a -> a -> Tree a -> Tree a
cut _ _ Nil = Nil
cut min max (Tree left val right)
  | val < min = cut min max right
  | val > max = cut min max left
  | otherwise = Tree (cut min max left) val (cut min max right)


testTree = Tree
  (Tree (Tree Nil 1 Nil) 3 (Tree (Tree Nil 4 Nil) 6 (Tree Nil 7 Nil)))
  8
  (Tree Nil 10 (Tree (Tree Nil 13 Nil) 14 Nil))
```

### Otočení v orientované sekvenci

Zdroj: [MFF Forum: Zkouška 12.7.2021](http://forum.matfyz.info/viewtopic.php?f=169&t=12193)

Na vstupu je daný seznam ``S`` obsahující dvojice ``(položka, orientace)``, kde
položky jsou obecné informace nějakého typu (například geny v chromozomu), a
orientace je typu ``Bool`` (pro sousměrně a protisměrně). Volání funkce
``otoceni S`` má vydat seznam všech výsledků ``[Vs]`` jako seznam seznamů
dvojic stejného typu, kde jeden výsledek vznikne otočením nějaké souvislé části
``S``, přičemž v otočené části změníte informaci o směru. Délka otočené části
je od ``1`` do délky ``S``, tj. otáčenou spojitou část vybíráte všemi možnými
způsoby.

1. Napište (obecný) typ funkce ``otoceni``
2. Napište funkci ``otoceni``
3. Pracovala by Vaše implementace funkce otoceni na nekonečném vstupním
   seznamu? Šla by napsat správná implementace pro nekonečný seznam? (Stačí
   myšlenka: proč ano nebo proč ne.)

Příklad:

```haskell
> otoceni [('a',True),('b',True),('c',False)]
[[('a',False),('b',True),('c',False)],[('a',True),('b',False),('c',False)],[('b',False),('a',False),('c',False)],[('a',True),('b',True),('c',True)],[('a',True),('c',True),('b',False)],[('c',True),('b',False),('a',False)]]
```

Řešení:

```haskell
split3 :: [a] -> [([a], [a], [a])]
split3 as =
  [ (xs, ys, zs)
  | (n, _) <- zip [0 ..] (as ++ [undefined])
  , let prefix = take n as
  , let zs     = drop n as
  , (k, _) <- zip [0 ..] (prefix ++ [undefined])
  , let xs = take k prefix
  , let ys = drop k prefix
  , not $ null ys
  ]

otoceni :: [(a, Bool)] -> [[(a, Bool)]]
otoceni ps = [ xs ++ map flipPair ys ++ zs | (xs, ys, zs) <- split3 ps ]
  where flipPair (x, bool) = (x, not bool)
```

### Převážení binárního stromu II

Je zadán binární strom s vnitřními vrcholy typu

```haskell
data Bt a = Void | Node (Bt a) a (Bt a)
```

Definujte funkci ``prevaz`` která projde strom a pro každý vnitřní vrchol
prohodí levý a pravý podstrom, pokud je ve vstupním stromě vlevo víc vrcholů
než vpravo.

Příklad:

```haskell
> prevaz (Node (Node (Node Void 'a' Void) 'b' Void) 'c' (Node Void 'ď Void))
Node (Node Void 'ď Void) 'c' (Node Void 'b' (Node Void 'a' Void))
```

1. Napište co nejobecnější typ funkce ``prevaz`` a použitých pomocných funkcí
2. Napište funkci ``prevaz.``
3. Využíváte někde volání lambda-funkce nebo funkce s neúplně zadanými argumenty?

Řešení:

```haskell
data Bt a = Void | Node (Bt a) a (Bt a)
  deriving (Eq, Show)

prevaz :: Bt a -> Bt a
prevaz = snd . prevaz'

prevaz' :: Bt a -> (Int, Bt a)
prevaz' Void = (0, Void)
prevaz' (Node left val right)
  | leftN > rightN = (count, Node rightAns val leftAns)
  | otherwise      = (count, Node leftAns val rightAns)
 where
  (leftN , leftAns ) = prevaz' left
  (rightN, rightAns) = prevaz' right
  count              = leftN + rightN + 1
```

### Diskvalifikováni sousedi

Dostanete vstupní graf ``G``, neorientovaný a bez ohodnocení. Vypusťte z něho
opakovaně všechny vrcholy, které mají méně sousedů než dané ``k``. Vydejte
zbylý graf a seznam vrcholů v poradí, jak jste je vypouštěli.

1. Definujte vhodný typ ``Graf`` a pro graf, který používáte v další definici,
   přičemž parametr ``a`` je označení vrcholů.
2. Definujte funkci ``centrumG :: Eq a => Graf a -> Int -> (Graf a, [a])`` pro
   požadovaný výpočet.

Řešení:

```haskell
data Graph a = Graph [(a, [a])]
  deriving (Eq, Show)

unfold :: (a -> Bool) -> (a -> (a, b)) -> a -> (a, [b])
unfold done step x
  | done x
  = (x, [])
  | otherwise
  = let (newX, newY ) = step x
        (ansX, ansYs) = unfold done step newX
    in  (ansX, newY : ansYs)

centrumG :: Eq a => Graph a -> Int -> (Graph a, [a])
centrumG (Graph ps) k = (Graph ansPs, bs)
 where
  findSmallDegree = filter (\p -> length (snd p) < k)
  degreeAtLeast ps = null (findSmallDegree ps)
  removeSmallDegree ps =
    let toRemove = fst $ head $ findSmallDegree ps
    in  ( [ (v, filter (/= toRemove) ns) | (v, ns) <- ps, v /= toRemove ]
        , toRemove
        )

  (ansPs, bs) = unfold degreeAtLeast removeSmallDegree ps
```


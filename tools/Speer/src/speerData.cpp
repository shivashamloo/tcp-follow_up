#include <ncbi_pch.hpp>
#include <iostream>
#include <math.h>
#include "speerData.hpp"



TAADataMap refFreqMap;
TPropertyMap speerPropertyMap;
TPairDistanceMap speerPairDistanceMap;

//  Data needed for relative-entropy scoring
void SpeerInitializeBackgroundFreqData() {

    /*  From Saikat's email...
    A  0.0780474700897585
    C  0.0192460110427568
    D  0.0536397361638076
    E  0.0629485981204668
    F  0.0385564048655621
    G  0.071551469592457    X wrong; see below for correct value
    H  0.0219922696262025
    I  0.0514196403000682
    K  0.0574383201866657
    L  0.090191394464413
    M  0.0224251883196316
    N  0.0448725775979007
    P  0.0520279465667327
    Q  0.0426458214465701
    R  0.0512953149316987
    S  0.0711984743501224
    T  0.0584129422708473
    V  0.0644094211988074
    W  0.013298374223799
    Y  0.0321647488738564
    */

    if (refFreqMap.size() == 0) {
        refFreqMap.insert(TAADataVT('a', 0.0780474700897585));
        refFreqMap.insert(TAADataVT('b', 0.0536397361638076 + 0.0448725775979007));  // sum of d, n
        refFreqMap.insert(TAADataVT('c', 0.0192460110427568));
        refFreqMap.insert(TAADataVT('d', 0.0536397361638076));
        refFreqMap.insert(TAADataVT('e', 0.0629485981204668));
        refFreqMap.insert(TAADataVT('f', 0.0385564048655621));
        refFreqMap.insert(TAADataVT('g', 0.0737715654561964));  
        refFreqMap.insert(TAADataVT('h', 0.0219922696262025));
        refFreqMap.insert(TAADataVT('i', 0.0514196403000682));
        refFreqMap.insert(TAADataVT('k', 0.0574383201866657));
        refFreqMap.insert(TAADataVT('l', 0.090191394464413));
        refFreqMap.insert(TAADataVT('m', 0.0224251883196316));
        refFreqMap.insert(TAADataVT('n', 0.0448725775979007));
        refFreqMap.insert(TAADataVT('p', 0.0520279465667327));
        refFreqMap.insert(TAADataVT('q', 0.0426458214465701));
        refFreqMap.insert(TAADataVT('r', 0.0512953149316987));
        refFreqMap.insert(TAADataVT('s', 0.0711984743501224));
        refFreqMap.insert(TAADataVT('t', 0.0584129422708473));
        refFreqMap.insert(TAADataVT('v', 0.0644094211988074));
        refFreqMap.insert(TAADataVT('w', 0.013298374223799));
        refFreqMap.insert(TAADataVT('x', 1.0));  //  any residue
        refFreqMap.insert(TAADataVT('y', 0.0321647488738564));
        refFreqMap.insert(TAADataVT('z', 0.0629485981204668 + 0.0426458214465701));  // sum of e, q
        refFreqMap.insert(TAADataVT('A', 0.0780474700897585));
        refFreqMap.insert(TAADataVT('B', 0.0536397361638076 + 0.0448725775979007));  // sum of D, N
        refFreqMap.insert(TAADataVT('C', 0.0192460110427568));
        refFreqMap.insert(TAADataVT('D', 0.0536397361638076));
        refFreqMap.insert(TAADataVT('E', 0.0629485981204668));
        refFreqMap.insert(TAADataVT('F', 0.0385564048655621));
        refFreqMap.insert(TAADataVT('G', 0.0737715654561964));  
        refFreqMap.insert(TAADataVT('H', 0.0219922696262025));
        refFreqMap.insert(TAADataVT('I', 0.0514196403000682));
        refFreqMap.insert(TAADataVT('K', 0.0574383201866657));
        refFreqMap.insert(TAADataVT('L', 0.090191394464413));
        refFreqMap.insert(TAADataVT('M', 0.0224251883196316));
        refFreqMap.insert(TAADataVT('N', 0.0448725775979007));
        refFreqMap.insert(TAADataVT('P', 0.0520279465667327));
        refFreqMap.insert(TAADataVT('Q', 0.0426458214465701));
        refFreqMap.insert(TAADataVT('R', 0.0512953149316987));
        refFreqMap.insert(TAADataVT('S', 0.0711984743501224));
        refFreqMap.insert(TAADataVT('T', 0.0584129422708473));
        refFreqMap.insert(TAADataVT('V', 0.0644094211988074));
        refFreqMap.insert(TAADataVT('W', 0.013298374223799));
        refFreqMap.insert(TAADataVT('X', 1.0));  //  any residue
        refFreqMap.insert(TAADataVT('Y', 0.0321647488738564));
        refFreqMap.insert(TAADataVT('Z', 0.0629485981204668 + 0.0426458214465701));  // sum of E, Q
        refFreqMap.insert(TAADataVT('-', 0.5));  //  gap
    }

/*
    for (TAADataIt it = refFreqMap.begin(); it != refFreqMap.end(); ++it) {
        cout << "refFreqMap[" << it->first << "] = " << it->second << endl;
    }
*/

}

//  Computes the score for the ambiguity "residues" B, Z and X using averages of upper-case residues.
//  Not using dataMap.insert(TAADataVT('x', sum)) to be sure to replace any existing value;
//  as such, this updates B, Z, X based on the current state of dataMap.
void SpeerComputeBZXPropertyValues(TAADataMap& dataMap) {
    unsigned int count = 0;
    double sum = 0.0;
    string::const_iterator cit, cend;
    TAADataMapCit mapEnd = dataMap.end();

    cend = AA_X.end();
    for ( cit = AA_X.begin(); cit != cend; ++cit) {
        if (dataMap.find(*cit) != mapEnd) {
            sum += dataMap[*cit];
        }
        ++count;  //  average over all residues, not just the ones defined w/ a value
    }
    if (count > 0) {
        sum /= count;
    } else {
        sum = 0.0;
    }
    dataMap['x'] = sum;  
    dataMap['X'] = sum;  

    count = 0;
    sum = 0.0;
    cend = AA_B.end();
    mapEnd = dataMap.end();
    for (cit = AA_B.begin(); cit != cend; ++cit) {
        if (dataMap.find(*cit) != mapEnd) {
            sum += dataMap[*cit];
        }
        ++count;  //  average over all residues, not just the ones defined w/ a value
    }
    if (count > 0) {
        sum /= count;
    } else {
        sum = 0.0;
    }
    dataMap['b'] = sum;  
    dataMap['B'] = sum;  


    count = 0;
    sum = 0.0;
    cend = AA_Z.end();
    mapEnd = dataMap.end();
    for (cit = AA_Z.begin(); cit != cend; ++cit) {
        if (dataMap.find(*cit) != mapEnd) {
            sum += dataMap[*cit];
        }
        ++count;  //  average over all residues, not just the ones defined w/ a value
    }
    if (count > 0) {
        sum /= count;
    } else {
        sum = 0.0;
    }
    dataMap['z'] = sum;  
    dataMap['Z'] = sum;  

}

void SpeerInitializeContactNumber() {

    /*  From ALL-INDEX
14a_contact_number
A 0.405737704918033 
R 0.332991803278688 
N 0.156762295081967 
D 0.00614754098360652 
C 0.905737704918033 
Q 0.145491803278689 
E 0.055327868852459 
G 0.262295081967213 
H 0.559426229508197 
I 1 
L 0.941598360655738 
K 0 
M 0.787909836065574 
F 0.968237704918033 
P 0.117827868852459 
S 0.137295081967213 
T 0.305327868852459 
W 0.961065573770492 
Y 0.648565573770492 
V 0.88422131147541 
    */

    TAADataMap& speerContactNumberMap = speerPropertyMap["14a_contact_number"];
    speerContactNumberMap.clear();

    speerContactNumberMap.insert(TAADataVT('a', 0.405737704918033 ));
    speerContactNumberMap.insert(TAADataVT('r', 0.332991803278688 ));
    speerContactNumberMap.insert(TAADataVT('n', 0.156762295081967 ));
    speerContactNumberMap.insert(TAADataVT('d', 0.00614754098360652 ));
    speerContactNumberMap.insert(TAADataVT('c', 0.905737704918033 ));
    speerContactNumberMap.insert(TAADataVT('q', 0.145491803278689 ));
    speerContactNumberMap.insert(TAADataVT('e', 0.055327868852459 ));
    speerContactNumberMap.insert(TAADataVT('g', 0.262295081967213 ));
    speerContactNumberMap.insert(TAADataVT('h', 0.559426229508197 ));
    speerContactNumberMap.insert(TAADataVT('i', 1 ));
    speerContactNumberMap.insert(TAADataVT('l', 0.941598360655738 ));
    speerContactNumberMap.insert(TAADataVT('k', 0 ));
    speerContactNumberMap.insert(TAADataVT('m', 0.787909836065574 ));
    speerContactNumberMap.insert(TAADataVT('f', 0.968237704918033 ));
    speerContactNumberMap.insert(TAADataVT('p', 0.117827868852459 ));
    speerContactNumberMap.insert(TAADataVT('s', 0.137295081967213 ));
    speerContactNumberMap.insert(TAADataVT('t', 0.305327868852459 ));
    speerContactNumberMap.insert(TAADataVT('w', 0.961065573770492 ));
    speerContactNumberMap.insert(TAADataVT('y', 0.648565573770492 ));
    speerContactNumberMap.insert(TAADataVT('v', 0.88422131147541 ));
    speerContactNumberMap.insert(TAADataVT('A', 0.405737704918033 ));
    speerContactNumberMap.insert(TAADataVT('R', 0.332991803278688 ));
    speerContactNumberMap.insert(TAADataVT('N', 0.156762295081967 ));
    speerContactNumberMap.insert(TAADataVT('D', 0.00614754098360652 ));
    speerContactNumberMap.insert(TAADataVT('C', 0.905737704918033 ));
    speerContactNumberMap.insert(TAADataVT('Q', 0.145491803278689 ));
    speerContactNumberMap.insert(TAADataVT('E', 0.055327868852459 ));
    speerContactNumberMap.insert(TAADataVT('G', 0.262295081967213 ));
    speerContactNumberMap.insert(TAADataVT('H', 0.559426229508197 ));
    speerContactNumberMap.insert(TAADataVT('I', 1 ));
    speerContactNumberMap.insert(TAADataVT('L', 0.941598360655738 ));
    speerContactNumberMap.insert(TAADataVT('K', 0 ));
    speerContactNumberMap.insert(TAADataVT('M', 0.787909836065574 ));
    speerContactNumberMap.insert(TAADataVT('F', 0.968237704918033 ));
    speerContactNumberMap.insert(TAADataVT('P', 0.117827868852459 ));
    speerContactNumberMap.insert(TAADataVT('S', 0.137295081967213 ));
    speerContactNumberMap.insert(TAADataVT('T', 0.305327868852459 ));
    speerContactNumberMap.insert(TAADataVT('W', 0.961065573770492 ));
    speerContactNumberMap.insert(TAADataVT('Y', 0.648565573770492 ));
    speerContactNumberMap.insert(TAADataVT('V', 0.88422131147541 ));
    speerContactNumberMap.insert(TAADataVT('-', -1));  //  gap

    SpeerComputeBZXPropertyValues(speerContactNumberMap);

/*
    for (TAADataIt it = speerContactNumberMap.begin(); it != speerContactNumberMap.end(); ++it) {
        cout << "speerContactNumberMap[" << it->first << "] = " << it->second << endl;
    }
*/

}

void SpeerInitializeAbsoluteEntropy() {

    /*  From ALL-INDEX
absolute_entropy
A 0.140535591668574 
R 1 
N 0.388189517051957 
D 0.36438544289311 
C 0.665827420462348 
Q 0.500801098649576 
E 0.463263904783703 
G 0 
H 0.944151979858091 
I 0.571526665140764 
L 0.592355230029755 
K 0.880521858548867 
M 0.699931334401465 
F 0.602426184481575 
P 0.331197070267796 
S 0.249713893339437 
T 0.269169146257725 
W 0.807049668116274 
Y 0.604486152437629 
V 0.412222476539254 
    */

    TAADataMap& speerAbsoluteEntropyMap = speerPropertyMap["absolute_entropy"];
    speerAbsoluteEntropyMap.clear();

    speerAbsoluteEntropyMap.insert(TAADataVT('a', 0.140535591668574 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('r', 1 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('n', 0.388189517051957 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('d', 0.36438544289311 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('c', 0.665827420462348 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('q', 0.500801098649576 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('e', 0.463263904783703 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('g', 0 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('h', 0.944151979858091 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('i', 0.571526665140764 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('l', 0.592355230029755 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('k', 0.880521858548867 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('m', 0.699931334401465 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('f', 0.602426184481575 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('p', 0.331197070267796 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('s', 0.249713893339437 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('t', 0.269169146257725 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('w', 0.807049668116274 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('y', 0.604486152437629 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('v', 0.412222476539254 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('A', 0.140535591668574 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('R', 1 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('N', 0.388189517051957 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('D', 0.36438544289311 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('C', 0.665827420462348 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('Q', 0.500801098649576 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('E', 0.463263904783703 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('G', 0 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('H', 0.944151979858091 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('I', 0.571526665140764 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('L', 0.592355230029755 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('K', 0.880521858548867 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('M', 0.699931334401465 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('F', 0.602426184481575 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('P', 0.331197070267796 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('S', 0.249713893339437 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('T', 0.269169146257725 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('W', 0.807049668116274 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('Y', 0.604486152437629 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('V', 0.412222476539254 ));
    speerAbsoluteEntropyMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerAbsoluteEntropyMap);

/*
  for (TAADataIt it = speerAbsoluteEntropyMap.begin(); it != speerAbsoluteEntropyMap.end(); ++it) {
  cout << "speerAbsoluteEntropyMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeAccSurfArea() {

    /*  From ALL-INDEX
accessible_surface_area
A 0.34499263622975 
R 0.921944035346097 
N 0.538659793814433 
D 0.525036818851252 
C 0.497790868924889 
Q 0.654270986745213 
E 0.673416789396171 
G 0.193667157584683 
H 0.692562592047128 
I 0.670839469808542 
L 0.639543446244477 
K 0.792341678939617 
M 0.727540500736377 
F 0.841678939617084 
P 0 
S 0.403166421207658 
T 0.52319587628866 
W 1 
Y 0.883284241531664 
V 0.578792341678939 
    */

    TAADataMap& speerAccSurfAreaMap = speerPropertyMap["accessible_surface_area"];
    speerAccSurfAreaMap.clear();

    speerAccSurfAreaMap.insert(TAADataVT('a', 0.34499263622975 ));
    speerAccSurfAreaMap.insert(TAADataVT('r', 0.921944035346097 ));
    speerAccSurfAreaMap.insert(TAADataVT('n', 0.538659793814433 ));
    speerAccSurfAreaMap.insert(TAADataVT('d', 0.525036818851252 ));
    speerAccSurfAreaMap.insert(TAADataVT('c', 0.497790868924889 ));
    speerAccSurfAreaMap.insert(TAADataVT('q', 0.654270986745213 ));
    speerAccSurfAreaMap.insert(TAADataVT('e', 0.673416789396171 ));
    speerAccSurfAreaMap.insert(TAADataVT('g', 0.193667157584683 ));
    speerAccSurfAreaMap.insert(TAADataVT('h', 0.692562592047128 ));
    speerAccSurfAreaMap.insert(TAADataVT('i', 0.670839469808542 ));
    speerAccSurfAreaMap.insert(TAADataVT('l', 0.639543446244477 ));
    speerAccSurfAreaMap.insert(TAADataVT('k', 0.792341678939617 ));
    speerAccSurfAreaMap.insert(TAADataVT('m', 0.727540500736377 ));
    speerAccSurfAreaMap.insert(TAADataVT('f', 0.841678939617084 ));
    speerAccSurfAreaMap.insert(TAADataVT('p', 0 ));
    speerAccSurfAreaMap.insert(TAADataVT('s', 0.403166421207658 ));
    speerAccSurfAreaMap.insert(TAADataVT('t', 0.52319587628866 ));
    speerAccSurfAreaMap.insert(TAADataVT('w', 1 ));
    speerAccSurfAreaMap.insert(TAADataVT('y', 0.883284241531664 ));
    speerAccSurfAreaMap.insert(TAADataVT('v', 0.578792341678939 ));
    speerAccSurfAreaMap.insert(TAADataVT('A', 0.34499263622975 ));
    speerAccSurfAreaMap.insert(TAADataVT('R', 0.921944035346097 ));
    speerAccSurfAreaMap.insert(TAADataVT('N', 0.538659793814433 ));
    speerAccSurfAreaMap.insert(TAADataVT('D', 0.525036818851252 ));
    speerAccSurfAreaMap.insert(TAADataVT('C', 0.497790868924889 ));
    speerAccSurfAreaMap.insert(TAADataVT('Q', 0.654270986745213 ));
    speerAccSurfAreaMap.insert(TAADataVT('E', 0.673416789396171 ));
    speerAccSurfAreaMap.insert(TAADataVT('G', 0.193667157584683 ));
    speerAccSurfAreaMap.insert(TAADataVT('H', 0.692562592047128 ));
    speerAccSurfAreaMap.insert(TAADataVT('I', 0.670839469808542 ));
    speerAccSurfAreaMap.insert(TAADataVT('L', 0.639543446244477 ));
    speerAccSurfAreaMap.insert(TAADataVT('K', 0.792341678939617 ));
    speerAccSurfAreaMap.insert(TAADataVT('M', 0.727540500736377 ));
    speerAccSurfAreaMap.insert(TAADataVT('F', 0.841678939617084 ));
    speerAccSurfAreaMap.insert(TAADataVT('P', 0 ));
    speerAccSurfAreaMap.insert(TAADataVT('S', 0.403166421207658 ));
    speerAccSurfAreaMap.insert(TAADataVT('T', 0.52319587628866 ));
    speerAccSurfAreaMap.insert(TAADataVT('W', 1 ));
    speerAccSurfAreaMap.insert(TAADataVT('Y', 0.883284241531664 ));
    speerAccSurfAreaMap.insert(TAADataVT('V', 0.578792341678939 ));
    speerAccSurfAreaMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerAccSurfAreaMap);

/*
  for (TAADataIt it = speerAccSurfAreaMap.begin(); it != speerAccSurfAreaMap.end(); ++it) {
  cout << "speerAccSurfAreaMap[" << it->first << "] = " << it->second << endl;
  }
*/

}
void SpeerInitializeHydrationPotl() {

    /*  From ALL-INDEX
        hydration_potential
A 0.97982967279247 
R 0 
N 0.458987001344689 
D 0.402061855670103 
C 0.837292693859256 
Q 0.472433886149709 
E 0.435679067682654 
G 1 
H 0.432541461228149 
I 0.989242492155984 
L 0.995069475571493 
K 0.466158673240699 
M 0.82653518601524 
F 0.858807709547288 
P 0.727924697445092 
S 0.666069027341999 
T 0.674137158225011 
W 0.629314208874944 
Y 0.619004930524429 
V 0.982070820259973 
    */

    TAADataMap& speerHydrationPotlMap = speerPropertyMap["hydration_potential"];
    speerHydrationPotlMap.clear();

    speerHydrationPotlMap.insert(TAADataVT('a', 0.97982967279247 ));
    speerHydrationPotlMap.insert(TAADataVT('r', 0 ));
    speerHydrationPotlMap.insert(TAADataVT('n', 0.458987001344689 ));
    speerHydrationPotlMap.insert(TAADataVT('d', 0.402061855670103 ));
    speerHydrationPotlMap.insert(TAADataVT('c', 0.837292693859256 ));
    speerHydrationPotlMap.insert(TAADataVT('q', 0.472433886149709 ));
    speerHydrationPotlMap.insert(TAADataVT('e', 0.435679067682654 ));
    speerHydrationPotlMap.insert(TAADataVT('g', 1 ));
    speerHydrationPotlMap.insert(TAADataVT('h', 0.432541461228149 ));
    speerHydrationPotlMap.insert(TAADataVT('i', 0.989242492155984 ));
    speerHydrationPotlMap.insert(TAADataVT('l', 0.995069475571493 ));
    speerHydrationPotlMap.insert(TAADataVT('k', 0.466158673240699 ));
    speerHydrationPotlMap.insert(TAADataVT('m', 0.82653518601524 ));
    speerHydrationPotlMap.insert(TAADataVT('f', 0.858807709547288 ));
    speerHydrationPotlMap.insert(TAADataVT('p', 0.727924697445092 ));
    speerHydrationPotlMap.insert(TAADataVT('s', 0.666069027341999 ));
    speerHydrationPotlMap.insert(TAADataVT('t', 0.674137158225011 ));
    speerHydrationPotlMap.insert(TAADataVT('w', 0.629314208874944 ));
    speerHydrationPotlMap.insert(TAADataVT('y', 0.619004930524429 ));
    speerHydrationPotlMap.insert(TAADataVT('v', 0.982070820259973 ));
    speerHydrationPotlMap.insert(TAADataVT('A', 0.97982967279247 ));
    speerHydrationPotlMap.insert(TAADataVT('R', 0 ));
    speerHydrationPotlMap.insert(TAADataVT('N', 0.458987001344689 ));
    speerHydrationPotlMap.insert(TAADataVT('D', 0.402061855670103 ));
    speerHydrationPotlMap.insert(TAADataVT('C', 0.837292693859256 ));
    speerHydrationPotlMap.insert(TAADataVT('Q', 0.472433886149709 ));
    speerHydrationPotlMap.insert(TAADataVT('E', 0.435679067682654 ));
    speerHydrationPotlMap.insert(TAADataVT('G', 1 ));
    speerHydrationPotlMap.insert(TAADataVT('H', 0.432541461228149 ));
    speerHydrationPotlMap.insert(TAADataVT('I', 0.989242492155984 ));
    speerHydrationPotlMap.insert(TAADataVT('L', 0.995069475571493 ));
    speerHydrationPotlMap.insert(TAADataVT('K', 0.466158673240699 ));
    speerHydrationPotlMap.insert(TAADataVT('M', 0.82653518601524 ));
    speerHydrationPotlMap.insert(TAADataVT('F', 0.858807709547288 ));
    speerHydrationPotlMap.insert(TAADataVT('P', 0.727924697445092 ));
    speerHydrationPotlMap.insert(TAADataVT('S', 0.666069027341999 ));
    speerHydrationPotlMap.insert(TAADataVT('T', 0.674137158225011 ));
    speerHydrationPotlMap.insert(TAADataVT('W', 0.629314208874944 ));
    speerHydrationPotlMap.insert(TAADataVT('Y', 0.619004930524429 ));
    speerHydrationPotlMap.insert(TAADataVT('V', 0.982070820259973 ));
    speerHydrationPotlMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerHydrationPotlMap);

/*
  for (TAADataIt it = speerHydrationPotlMap.begin(); it != speerHydrationPotlMap.end(); ++it) {
  cout << "speerHydrationPotlMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeHydropathyIndex() {

    /*  From ALL-INDEX
hydropathy_index
A 0.7 
R 0 
N 0.111111111111111 
D 0.111111111111111 
C 0.777777777777778 
Q 0.111111111111111 
E 0.111111111111111 
G 0.455555555555555 
H 0.144444444444444 
I 1 
L 0.922222222222222 
K 0.0666666666666667 
M 0.711111111111111 
F 0.811111111111111 
P 0.322222222222222 
S 0.411111111111111 
T 0.422222222222222 
W 0.4 
Y 0.355555555555556 
V 0.966666666666667 
    */

    TAADataMap& speerHydropathyIndexMap = speerPropertyMap["hydropathy_index"];
    speerHydropathyIndexMap.clear();

    speerHydropathyIndexMap.insert(TAADataVT('a', 0.7 ));
    speerHydropathyIndexMap.insert(TAADataVT('r', 0 ));
    speerHydropathyIndexMap.insert(TAADataVT('n', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('d', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('c', 0.777777777777778 ));
    speerHydropathyIndexMap.insert(TAADataVT('q', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('e', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('g', 0.455555555555555 ));
    speerHydropathyIndexMap.insert(TAADataVT('h', 0.144444444444444 ));
    speerHydropathyIndexMap.insert(TAADataVT('i', 1 ));
    speerHydropathyIndexMap.insert(TAADataVT('l', 0.922222222222222 ));
    speerHydropathyIndexMap.insert(TAADataVT('k', 0.0666666666666667 ));
    speerHydropathyIndexMap.insert(TAADataVT('m', 0.711111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('f', 0.811111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('p', 0.322222222222222 ));
    speerHydropathyIndexMap.insert(TAADataVT('s', 0.411111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('t', 0.422222222222222 ));
    speerHydropathyIndexMap.insert(TAADataVT('w', 0.4 ));
    speerHydropathyIndexMap.insert(TAADataVT('y', 0.355555555555556 ));
    speerHydropathyIndexMap.insert(TAADataVT('v', 0.966666666666667 ));
    speerHydropathyIndexMap.insert(TAADataVT('A', 0.7 ));
    speerHydropathyIndexMap.insert(TAADataVT('R', 0 ));
    speerHydropathyIndexMap.insert(TAADataVT('N', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('D', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('C', 0.777777777777778 ));
    speerHydropathyIndexMap.insert(TAADataVT('Q', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('E', 0.111111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('G', 0.455555555555555 ));
    speerHydropathyIndexMap.insert(TAADataVT('H', 0.144444444444444 ));
    speerHydropathyIndexMap.insert(TAADataVT('I', 1 ));
    speerHydropathyIndexMap.insert(TAADataVT('L', 0.922222222222222 ));
    speerHydropathyIndexMap.insert(TAADataVT('K', 0.0666666666666667 ));
    speerHydropathyIndexMap.insert(TAADataVT('M', 0.711111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('F', 0.811111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('P', 0.322222222222222 ));
    speerHydropathyIndexMap.insert(TAADataVT('S', 0.411111111111111 ));
    speerHydropathyIndexMap.insert(TAADataVT('T', 0.422222222222222 ));
    speerHydropathyIndexMap.insert(TAADataVT('W', 0.4 ));
    speerHydropathyIndexMap.insert(TAADataVT('Y', 0.355555555555556 ));
    speerHydropathyIndexMap.insert(TAADataVT('V', 0.966666666666667 ));
    speerHydropathyIndexMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerHydropathyIndexMap);

/*
  for (TAADataIt it = speerHydropathyIndexMap.begin(); it != speerHydropathyIndexMap.end(); ++it) {
  cout << "speerHydropathyIndexMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeHydrophobicityIndex() {

    /*  From ALL-INDEX
        hydrophobicity_index
A 0.230188679245283 
R 0.226415094339623 
N 0.0226415094339623 
D 0.173584905660377 
C 0.40377358490566 
Q 0 
E 0.177358490566038 
G 0.0264150943396226 
H 0.230188679245283 
I 0.837735849056604 
L 0.577358490566038 
K 0.433962264150943 
M 0.445283018867925 
F 0.762264150943396 
P 0.735849056603774 
S 0.0188679245283019 
T 0.0188679245283019 
W 1 
Y 0.709433962264151 
V 0.49811320754717 
    */

    TAADataMap& speerHydrophobicityIndexMap = speerPropertyMap["hydrophobicity_index"];
    speerHydrophobicityIndexMap.clear();

    speerHydrophobicityIndexMap.insert(TAADataVT('a', 0.230188679245283 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('r', 0.226415094339623 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('n', 0.0226415094339623 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('d', 0.173584905660377 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('c', 0.40377358490566 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('q', 0 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('e', 0.177358490566038 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('g', 0.0264150943396226 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('h', 0.230188679245283 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('i', 0.837735849056604 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('l', 0.577358490566038 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('k', 0.433962264150943 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('m', 0.445283018867925 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('f', 0.762264150943396 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('p', 0.735849056603774 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('s', 0.0188679245283019 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('t', 0.0188679245283019 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('w', 1 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('y', 0.709433962264151 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('v', 0.49811320754717 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('A', 0.230188679245283 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('R', 0.226415094339623 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('N', 0.0226415094339623 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('D', 0.173584905660377 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('C', 0.40377358490566 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('Q', 0 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('E', 0.177358490566038 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('G', 0.0264150943396226 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('H', 0.230188679245283 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('I', 0.837735849056604 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('L', 0.577358490566038 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('K', 0.433962264150943 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('M', 0.445283018867925 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('F', 0.762264150943396 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('P', 0.735849056603774 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('S', 0.0188679245283019 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('T', 0.0188679245283019 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('W', 1 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('Y', 0.709433962264151 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('V', 0.49811320754717 ));
    speerHydrophobicityIndexMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerHydrophobicityIndexMap);

/*
  for (TAADataIt it = speerHydrophobicityIndexMap.begin(); it != speerHydrophobicityIndexMap.end(); ++it) {
  cout << "speerHydrophobicityIndexMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeIsoElecPoint() {

    /*  From ALL-INDEX
isoelectric_point
A 0.404255319148936 
R 1 
N 0.330413016270338 
D 0 
C 0.285356695869837 
Q 0.360450563204005 
E 0.0563204005006258 
G 0.400500625782228 
H 0.603254067584481 
I 0.406758448060075 
L 0.401752190237797 
K 0.872340425531915 
M 0.37171464330413 
F 0.339173967459324 
P 0.44180225281602 
S 0.364205256570713 
T 0.361702127659574 
W 0.390488110137672 
Y 0.361702127659574 
V 0.399249061326658 
    */

    TAADataMap& speerIsoElecPointMap = speerPropertyMap["isoelectric_point"];
    speerIsoElecPointMap.clear();

    speerIsoElecPointMap.insert(TAADataVT('a', 0.404255319148936 ));
    speerIsoElecPointMap.insert(TAADataVT('r', 1 ));
    speerIsoElecPointMap.insert(TAADataVT('n', 0.330413016270338 ));
    speerIsoElecPointMap.insert(TAADataVT('d', 0 ));
    speerIsoElecPointMap.insert(TAADataVT('c', 0.285356695869837 ));
    speerIsoElecPointMap.insert(TAADataVT('q', 0.360450563204005 ));
    speerIsoElecPointMap.insert(TAADataVT('e', 0.0563204005006258 ));
    speerIsoElecPointMap.insert(TAADataVT('g', 0.400500625782228 ));
    speerIsoElecPointMap.insert(TAADataVT('h', 0.603254067584481 ));
    speerIsoElecPointMap.insert(TAADataVT('i', 0.406758448060075 ));
    speerIsoElecPointMap.insert(TAADataVT('l', 0.401752190237797 ));
    speerIsoElecPointMap.insert(TAADataVT('k', 0.872340425531915 ));
    speerIsoElecPointMap.insert(TAADataVT('m', 0.37171464330413 ));
    speerIsoElecPointMap.insert(TAADataVT('f', 0.339173967459324 ));
    speerIsoElecPointMap.insert(TAADataVT('p', 0.44180225281602 ));
    speerIsoElecPointMap.insert(TAADataVT('s', 0.364205256570713 ));
    speerIsoElecPointMap.insert(TAADataVT('t', 0.361702127659574 ));
    speerIsoElecPointMap.insert(TAADataVT('w', 0.390488110137672 ));
    speerIsoElecPointMap.insert(TAADataVT('y', 0.361702127659574 ));
    speerIsoElecPointMap.insert(TAADataVT('v', 0.399249061326658 ));
    speerIsoElecPointMap.insert(TAADataVT('A', 0.404255319148936 ));
    speerIsoElecPointMap.insert(TAADataVT('R', 1 ));
    speerIsoElecPointMap.insert(TAADataVT('N', 0.330413016270338 ));
    speerIsoElecPointMap.insert(TAADataVT('D', 0 ));
    speerIsoElecPointMap.insert(TAADataVT('C', 0.285356695869837 ));
    speerIsoElecPointMap.insert(TAADataVT('Q', 0.360450563204005 ));
    speerIsoElecPointMap.insert(TAADataVT('E', 0.0563204005006258 ));
    speerIsoElecPointMap.insert(TAADataVT('G', 0.400500625782228 ));
    speerIsoElecPointMap.insert(TAADataVT('H', 0.603254067584481 ));
    speerIsoElecPointMap.insert(TAADataVT('I', 0.406758448060075 ));
    speerIsoElecPointMap.insert(TAADataVT('L', 0.401752190237797 ));
    speerIsoElecPointMap.insert(TAADataVT('K', 0.872340425531915 ));
    speerIsoElecPointMap.insert(TAADataVT('M', 0.37171464330413 ));
    speerIsoElecPointMap.insert(TAADataVT('F', 0.339173967459324 ));
    speerIsoElecPointMap.insert(TAADataVT('P', 0.44180225281602 ));
    speerIsoElecPointMap.insert(TAADataVT('S', 0.364205256570713 ));
    speerIsoElecPointMap.insert(TAADataVT('T', 0.361702127659574 ));
    speerIsoElecPointMap.insert(TAADataVT('W', 0.390488110137672 ));
    speerIsoElecPointMap.insert(TAADataVT('Y', 0.361702127659574 ));
    speerIsoElecPointMap.insert(TAADataVT('V', 0.399249061326658 ));
    speerIsoElecPointMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerIsoElecPointMap);

/*
  for (TAADataIt it = speerIsoElecPointMap.begin(); it != speerIsoElecPointMap.end(); ++it) {
  cout << "speerIsoElecPointMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeMolWeight() {

    /*  From ALL-INDEX
molecular_weight
A 0.108539134473949 
R 0.767438259657815 
N 0.441666021522025 
D 0.449252922505226 
C 0.356739180924363 
Q 0.550282573352946 
E 0.557869474336146 
G 0 
H 0.620035611984207 
I 0.434311372609739 
L 0.434311372609739 
K 0.550592242780831 
M 0.573972284586204 
F 0.697685221026554 
P 0.310133932027561 
S 0.232406905628242 
T 0.341023457459162 
W 1 
Y 0.821552992180847 
V 0.32577223813579 
    */

    TAADataMap& speerMolWeightMap = speerPropertyMap["molecular_weight"];
    speerMolWeightMap.clear();

    speerMolWeightMap.insert(TAADataVT('a', 0.108539134473949 ));
    speerMolWeightMap.insert(TAADataVT('r', 0.767438259657815 ));
    speerMolWeightMap.insert(TAADataVT('n', 0.441666021522025 ));
    speerMolWeightMap.insert(TAADataVT('d', 0.449252922505226 ));
    speerMolWeightMap.insert(TAADataVT('c', 0.356739180924363 ));
    speerMolWeightMap.insert(TAADataVT('q', 0.550282573352946 ));
    speerMolWeightMap.insert(TAADataVT('e', 0.557869474336146 ));
    speerMolWeightMap.insert(TAADataVT('g', 0 ));
    speerMolWeightMap.insert(TAADataVT('h', 0.620035611984207 ));
    speerMolWeightMap.insert(TAADataVT('i', 0.434311372609739 ));
    speerMolWeightMap.insert(TAADataVT('l', 0.434311372609739 ));
    speerMolWeightMap.insert(TAADataVT('k', 0.550592242780831 ));
    speerMolWeightMap.insert(TAADataVT('m', 0.573972284586204 ));
    speerMolWeightMap.insert(TAADataVT('f', 0.697685221026554 ));
    speerMolWeightMap.insert(TAADataVT('p', 0.310133932027561 ));
    speerMolWeightMap.insert(TAADataVT('s', 0.232406905628242 ));
    speerMolWeightMap.insert(TAADataVT('t', 0.341023457459162 ));
    speerMolWeightMap.insert(TAADataVT('w', 1 ));
    speerMolWeightMap.insert(TAADataVT('y', 0.821552992180847 ));
    speerMolWeightMap.insert(TAADataVT('v', 0.32577223813579 ));
    speerMolWeightMap.insert(TAADataVT('A', 0.108539134473949 ));
    speerMolWeightMap.insert(TAADataVT('R', 0.767438259657815 ));
    speerMolWeightMap.insert(TAADataVT('N', 0.441666021522025 ));
    speerMolWeightMap.insert(TAADataVT('D', 0.449252922505226 ));
    speerMolWeightMap.insert(TAADataVT('C', 0.356739180924363 ));
    speerMolWeightMap.insert(TAADataVT('Q', 0.550282573352946 ));
    speerMolWeightMap.insert(TAADataVT('E', 0.557869474336146 ));
    speerMolWeightMap.insert(TAADataVT('G', 0 ));
    speerMolWeightMap.insert(TAADataVT('H', 0.620035611984207 ));
    speerMolWeightMap.insert(TAADataVT('I', 0.434311372609739 ));
    speerMolWeightMap.insert(TAADataVT('L', 0.434311372609739 ));
    speerMolWeightMap.insert(TAADataVT('K', 0.550592242780831 ));
    speerMolWeightMap.insert(TAADataVT('M', 0.573972284586204 ));
    speerMolWeightMap.insert(TAADataVT('F', 0.697685221026554 ));
    speerMolWeightMap.insert(TAADataVT('P', 0.310133932027561 ));
    speerMolWeightMap.insert(TAADataVT('S', 0.232406905628242 ));
    speerMolWeightMap.insert(TAADataVT('T', 0.341023457459162 ));
    speerMolWeightMap.insert(TAADataVT('W', 1 ));
    speerMolWeightMap.insert(TAADataVT('Y', 0.821552992180847 ));
    speerMolWeightMap.insert(TAADataVT('V', 0.32577223813579 ));
    speerMolWeightMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerMolWeightMap);

/*
  for (TAADataIt it = speerMolWeightMap.begin(); it != speerMolWeightMap.end(); ++it) {
  cout << "speerMolWeightMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeNetCharge() {

    /*  From ALL-INDEX
net_charge
A 0.5 
R 1 
N 0.5 
D 0 
C 0.5 
Q 0.5 
E 0 
G 0.5 
H 0.5 
I 0.5 
L 0.5 
K 1 
M 0.5 
F 0.5 
P 0.5 
S 0.5 
T 0.5 
W 0.5 
Y 0.5 
V 0.5 
    */

    TAADataMap& speerNetChargeMap = speerPropertyMap["net_charge"];
    speerNetChargeMap.clear();

    speerNetChargeMap.insert(TAADataVT('a', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('r', 1 ));
    speerNetChargeMap.insert(TAADataVT('n', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('d', 0 ));
    speerNetChargeMap.insert(TAADataVT('c', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('q', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('e', 0 ));
    speerNetChargeMap.insert(TAADataVT('g', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('h', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('i', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('l', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('k', 1 ));
    speerNetChargeMap.insert(TAADataVT('m', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('f', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('p', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('s', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('t', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('w', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('y', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('v', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('A', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('R', 1 ));
    speerNetChargeMap.insert(TAADataVT('N', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('D', 0 ));
    speerNetChargeMap.insert(TAADataVT('C', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('Q', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('E', 0 ));
    speerNetChargeMap.insert(TAADataVT('G', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('H', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('I', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('L', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('K', 1 ));
    speerNetChargeMap.insert(TAADataVT('M', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('F', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('P', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('S', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('T', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('W', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('Y', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('V', 0.5 ));
    speerNetChargeMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerNetChargeMap);

/*
  for (TAADataIt it = speerNetChargeMap.begin(); it != speerNetChargeMap.end(); ++it) {
  cout << "speerNetChargeMap[" << it->first << "] = " << it->second << endl;
  }
*/

}
void SpeerInitializeNormFlexParam() {

    /*  From ALL-INDEX
normalized_flexibility_parameters
A 0.404040404040404 
R 0.525252525252525 
N 0.727272727272727 
D 0.828282828282828 
C 0.0101010101010101 
Q 0.671717171717171 
E 0.95959595959596 
G 0.641414141414141 
H 0.232323232323232 
I 0.116161616161616 
L 0.156565656565657 
K 1 
M 0.242424242424242 
F 0.0555555555555556 
P 0.732323232323232 
S 0.717171717171717 
T 0.469696969696969 
W 0 
Y 0.126262626262626 
V 0.136363636363636 
    */

    TAADataMap& speerNormFlexParamMap = speerPropertyMap["normalized_flexibility_parameters"];
    speerNormFlexParamMap.clear();

    speerNormFlexParamMap.insert(TAADataVT('a', 0.404040404040404 ));
    speerNormFlexParamMap.insert(TAADataVT('r', 0.525252525252525 ));
    speerNormFlexParamMap.insert(TAADataVT('n', 0.727272727272727 ));
    speerNormFlexParamMap.insert(TAADataVT('d', 0.828282828282828 ));
    speerNormFlexParamMap.insert(TAADataVT('c', 0.0101010101010101 ));
    speerNormFlexParamMap.insert(TAADataVT('q', 0.671717171717171 ));
    speerNormFlexParamMap.insert(TAADataVT('e', 0.95959595959596 ));
    speerNormFlexParamMap.insert(TAADataVT('g', 0.641414141414141 ));
    speerNormFlexParamMap.insert(TAADataVT('h', 0.232323232323232 ));
    speerNormFlexParamMap.insert(TAADataVT('i', 0.116161616161616 ));
    speerNormFlexParamMap.insert(TAADataVT('l', 0.156565656565657 ));
    speerNormFlexParamMap.insert(TAADataVT('k', 1 ));
    speerNormFlexParamMap.insert(TAADataVT('m', 0.242424242424242 ));
    speerNormFlexParamMap.insert(TAADataVT('f', 0.0555555555555556 ));
    speerNormFlexParamMap.insert(TAADataVT('p', 0.732323232323232 ));
    speerNormFlexParamMap.insert(TAADataVT('s', 0.717171717171717 ));
    speerNormFlexParamMap.insert(TAADataVT('t', 0.469696969696969 ));
    speerNormFlexParamMap.insert(TAADataVT('w', 0 ));
    speerNormFlexParamMap.insert(TAADataVT('y', 0.126262626262626 ));
    speerNormFlexParamMap.insert(TAADataVT('v', 0.136363636363636 ));
    speerNormFlexParamMap.insert(TAADataVT('A', 0.404040404040404 ));
    speerNormFlexParamMap.insert(TAADataVT('R', 0.525252525252525 ));
    speerNormFlexParamMap.insert(TAADataVT('N', 0.727272727272727 ));
    speerNormFlexParamMap.insert(TAADataVT('D', 0.828282828282828 ));
    speerNormFlexParamMap.insert(TAADataVT('C', 0.0101010101010101 ));
    speerNormFlexParamMap.insert(TAADataVT('Q', 0.671717171717171 ));
    speerNormFlexParamMap.insert(TAADataVT('E', 0.95959595959596 ));
    speerNormFlexParamMap.insert(TAADataVT('G', 0.641414141414141 ));
    speerNormFlexParamMap.insert(TAADataVT('H', 0.232323232323232 ));
    speerNormFlexParamMap.insert(TAADataVT('I', 0.116161616161616 ));
    speerNormFlexParamMap.insert(TAADataVT('L', 0.156565656565657 ));
    speerNormFlexParamMap.insert(TAADataVT('K', 1 ));
    speerNormFlexParamMap.insert(TAADataVT('M', 0.242424242424242 ));
    speerNormFlexParamMap.insert(TAADataVT('F', 0.0555555555555556 ));
    speerNormFlexParamMap.insert(TAADataVT('P', 0.732323232323232 ));
    speerNormFlexParamMap.insert(TAADataVT('S', 0.717171717171717 ));
    speerNormFlexParamMap.insert(TAADataVT('T', 0.469696969696969 ));
    speerNormFlexParamMap.insert(TAADataVT('W', 0 ));
    speerNormFlexParamMap.insert(TAADataVT('Y', 0.126262626262626 ));
    speerNormFlexParamMap.insert(TAADataVT('V', 0.136363636363636 ));
    speerNormFlexParamMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerNormFlexParamMap);

/*
    for (TAADataIt it = speerNormFlexParamMap.begin(); it != speerNormFlexParamMap.end(); ++it) {
        cout << "speerNormFlexParamMap[" << it->first << "] = " << it->second << endl;
    }
*/

}

void SpeerInitializeRelativeMutability() {

    /*  From ALL-INDEX
relative_mutability
A 0.815217391304348 
R 0.630434782608696 
N 0.858695652173913 
D 0.66304347826087 
C 0.206521739130435 
Q 0.641304347826087 
E 0.565217391304348 
G 0.271739130434783 
H 0.717391304347826 
I 0.847826086956522 
L 0.315217391304348 
K 0.510869565217391 
M 0.739130434782609 
F 0.282608695652174 
P 0.358695652173913 
S 1 
T 0.891304347826087 
W 0 
Y 0.271739130434783 
V 0.793478260869565 
    */

    TAADataMap& speerRelativeMutabilityMap = speerPropertyMap["relative_mutability"];
    speerRelativeMutabilityMap.clear();

    speerRelativeMutabilityMap.insert(TAADataVT('a', 0.815217391304348 ));
    speerRelativeMutabilityMap.insert(TAADataVT('r', 0.630434782608696 ));
    speerRelativeMutabilityMap.insert(TAADataVT('n', 0.858695652173913 ));
    speerRelativeMutabilityMap.insert(TAADataVT('d', 0.66304347826087 ));
    speerRelativeMutabilityMap.insert(TAADataVT('c', 0.206521739130435 ));
    speerRelativeMutabilityMap.insert(TAADataVT('q', 0.641304347826087 ));
    speerRelativeMutabilityMap.insert(TAADataVT('e', 0.565217391304348 ));
    speerRelativeMutabilityMap.insert(TAADataVT('g', 0.271739130434783 ));
    speerRelativeMutabilityMap.insert(TAADataVT('h', 0.717391304347826 ));
    speerRelativeMutabilityMap.insert(TAADataVT('i', 0.847826086956522 ));
    speerRelativeMutabilityMap.insert(TAADataVT('l', 0.315217391304348 ));
    speerRelativeMutabilityMap.insert(TAADataVT('k', 0.510869565217391 ));
    speerRelativeMutabilityMap.insert(TAADataVT('m', 0.739130434782609 ));
    speerRelativeMutabilityMap.insert(TAADataVT('f', 0.282608695652174 ));
    speerRelativeMutabilityMap.insert(TAADataVT('p', 0.358695652173913 ));
    speerRelativeMutabilityMap.insert(TAADataVT('s', 1 ));
    speerRelativeMutabilityMap.insert(TAADataVT('t', 0.891304347826087 ));
    speerRelativeMutabilityMap.insert(TAADataVT('w', 0 ));
    speerRelativeMutabilityMap.insert(TAADataVT('y', 0.271739130434783 ));
    speerRelativeMutabilityMap.insert(TAADataVT('v', 0.793478260869565 ));
    speerRelativeMutabilityMap.insert(TAADataVT('A', 0.815217391304348 ));
    speerRelativeMutabilityMap.insert(TAADataVT('R', 0.630434782608696 ));
    speerRelativeMutabilityMap.insert(TAADataVT('N', 0.858695652173913 ));
    speerRelativeMutabilityMap.insert(TAADataVT('D', 0.66304347826087 ));
    speerRelativeMutabilityMap.insert(TAADataVT('C', 0.206521739130435 ));
    speerRelativeMutabilityMap.insert(TAADataVT('Q', 0.641304347826087 ));
    speerRelativeMutabilityMap.insert(TAADataVT('E', 0.565217391304348 ));
    speerRelativeMutabilityMap.insert(TAADataVT('G', 0.271739130434783 ));
    speerRelativeMutabilityMap.insert(TAADataVT('H', 0.717391304347826 ));
    speerRelativeMutabilityMap.insert(TAADataVT('I', 0.847826086956522 ));
    speerRelativeMutabilityMap.insert(TAADataVT('L', 0.315217391304348 ));
    speerRelativeMutabilityMap.insert(TAADataVT('K', 0.510869565217391 ));
    speerRelativeMutabilityMap.insert(TAADataVT('M', 0.739130434782609 ));
    speerRelativeMutabilityMap.insert(TAADataVT('F', 0.282608695652174 ));
    speerRelativeMutabilityMap.insert(TAADataVT('P', 0.358695652173913 ));
    speerRelativeMutabilityMap.insert(TAADataVT('S', 1 ));
    speerRelativeMutabilityMap.insert(TAADataVT('T', 0.891304347826087 ));
    speerRelativeMutabilityMap.insert(TAADataVT('W', 0 ));
    speerRelativeMutabilityMap.insert(TAADataVT('Y', 0.271739130434783 ));
    speerRelativeMutabilityMap.insert(TAADataVT('V', 0.793478260869565 ));
    speerRelativeMutabilityMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerRelativeMutabilityMap);

/*
    for (TAADataIt it = speerRelativeMutabilityMap.begin(); it != speerRelativeMutabilityMap.end(); ++it) {
        cout << "speerRelativeMutabilityMap[" << it->first << "] = " << it->second << endl;
    }
*/

}
void SpeerInitializeAccSurfAreaInFoldedProtein() {

    /*  From ALL-INDEX
residue_accessible_surface_area_in_folded_protein
A 0.0886075949367089 
R 0.911392405063291 
N 0.569620253164557 
D 0.405063291139241 
C 0.0126582278481013 
Q 0.670886075949367 
E 0.392405063291139 
G 0.0632911392405063 
H 0.316455696202532 
I 0 
L 0.0632911392405063 
K 1 
M 0.164556962025316 
F 0.0759493670886076 
P 0.405063291139241 
S 0.329113924050633 
T 0.367088607594937 
W 0.177215189873418 
Y 0.531645569620253 
V 0 
    */

    TAADataMap& speerAccSurfAreaInFoldedProteinMap = speerPropertyMap["residue_accessible_surface_area_in_folded_protein"];
    speerAccSurfAreaInFoldedProteinMap.clear();

    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('a', 0.0886075949367089 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('r', 0.911392405063291 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('n', 0.569620253164557 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('d', 0.405063291139241 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('c', 0.0126582278481013 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('q', 0.670886075949367 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('e', 0.392405063291139 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('g', 0.0632911392405063 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('h', 0.316455696202532 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('i', 0 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('l', 0.0632911392405063 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('k', 1 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('m', 0.164556962025316 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('f', 0.0759493670886076 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('p', 0.405063291139241 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('s', 0.329113924050633 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('t', 0.367088607594937 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('w', 0.177215189873418 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('y', 0.531645569620253 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('v', 0 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('A', 0.0886075949367089 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('R', 0.911392405063291 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('N', 0.569620253164557 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('D', 0.405063291139241 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('C', 0.0126582278481013 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('Q', 0.670886075949367 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('E', 0.392405063291139 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('G', 0.0632911392405063 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('H', 0.316455696202532 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('I', 0 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('L', 0.0632911392405063 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('K', 1 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('M', 0.164556962025316 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('F', 0.0759493670886076 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('P', 0.405063291139241 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('S', 0.329113924050633 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('T', 0.367088607594937 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('W', 0.177215189873418 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('Y', 0.531645569620253 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('V', 0 ));
    speerAccSurfAreaInFoldedProteinMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerAccSurfAreaInFoldedProteinMap);

/*
  for (TAADataIt it = speerAccSurfAreaInFoldedProteinMap.begin(); it != speerAccSurfAreaInFoldedProteinMap.end(); ++it) {
  cout << "speerAccSurfAreaInFoldedProteinMap[" << it->first << "] = " << it->second << endl;
  }
*/

}

void SpeerInitializeSideChainOrientPref() {

    /*  From ALL-INDEX
side_chain_orientational_preference
A 0.217142857142857 
R 0.377142857142857 
N 0.448571428571429 
D 0.645714285714286 
C 0.0285714285714286 
Q 0.991428571428571 
E 0.571428571428571 
G 0.351428571428571 
H 0.131428571428571 
I 0.0485714285714286 
L 0.0314285714285714 
K 1 
M 0 
F 0.00857142857142856 
P 0.468571428571429 
S 0.345714285714286 
T 0.308571428571429 
W 0.1 
Y 0.377142857142857 
V 0.0542857142857143 
    */

    TAADataMap& speerSideChainOrientPrefMap = speerPropertyMap["side_chain_orientational_preference"];
    speerSideChainOrientPrefMap.clear();

    speerSideChainOrientPrefMap.insert(TAADataVT('a', 0.217142857142857 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('r', 0.377142857142857 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('n', 0.448571428571429 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('d', 0.645714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('c', 0.0285714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('q', 0.991428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('e', 0.571428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('g', 0.351428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('h', 0.131428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('i', 0.0485714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('l', 0.0314285714285714 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('k', 1 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('m', 0 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('f', 0.00857142857142856 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('p', 0.468571428571429 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('s', 0.345714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('t', 0.308571428571429 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('w', 0.1 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('y', 0.377142857142857 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('v', 0.0542857142857143 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('A', 0.217142857142857 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('R', 0.377142857142857 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('N', 0.448571428571429 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('D', 0.645714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('C', 0.0285714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('Q', 0.991428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('E', 0.571428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('G', 0.351428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('H', 0.131428571428571 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('I', 0.0485714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('L', 0.0314285714285714 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('K', 1 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('M', 0 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('F', 0.00857142857142856 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('P', 0.468571428571429 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('S', 0.345714285714286 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('T', 0.308571428571429 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('W', 0.1 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('Y', 0.377142857142857 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('V', 0.0542857142857143 ));
    speerSideChainOrientPrefMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerSideChainOrientPrefMap);

/*
    for (TAADataIt it = speerSideChainOrientPrefMap.begin(); it != speerSideChainOrientPrefMap.end(); ++it) {
        cout << "speerSideChainOrientPrefMap[" << it->first << "] = " << it->second << endl;
    }
*/

}

void SpeerInitializeOccurenceFreq() {

    /*  From ALL-INDEX
occurence_frequency
A 0.818181818181818 
R 0.480519480519481 
N 0.376623376623377 
D 0.493506493506494 
C 0.0779220779220779 
Q 0.350649350649351 
E 0.623376623376623 
G 0.779220779220779 
H 0.116883116883117 
I 0.506493506493506 
L 1 
K 0.584415584415584 
M 0.12987012987013 
F 0.337662337662338 
P 0.480519480519481 
S 0.714285714285714 
T 0.584415584415584 
W 0 
Y 0.233766233766234 
V 0.675324675324675 
    */

    TAADataMap& speerOccurenceFreqMap = speerPropertyMap["occurence_frequency"];
    speerOccurenceFreqMap.clear();

    speerOccurenceFreqMap.insert(TAADataVT('a', 0.818181818181818 ));
    speerOccurenceFreqMap.insert(TAADataVT('r', 0.480519480519481 ));
    speerOccurenceFreqMap.insert(TAADataVT('n', 0.376623376623377 ));
    speerOccurenceFreqMap.insert(TAADataVT('d', 0.493506493506494 ));
    speerOccurenceFreqMap.insert(TAADataVT('c', 0.0779220779220779 ));
    speerOccurenceFreqMap.insert(TAADataVT('q', 0.350649350649351 ));
    speerOccurenceFreqMap.insert(TAADataVT('e', 0.623376623376623 ));
    speerOccurenceFreqMap.insert(TAADataVT('g', 0.779220779220779 ));
    speerOccurenceFreqMap.insert(TAADataVT('h', 0.116883116883117 ));
    speerOccurenceFreqMap.insert(TAADataVT('i', 0.506493506493506 ));
    speerOccurenceFreqMap.insert(TAADataVT('l', 1 ));
    speerOccurenceFreqMap.insert(TAADataVT('k', 0.584415584415584 ));
    speerOccurenceFreqMap.insert(TAADataVT('m', 0.12987012987013 ));
    speerOccurenceFreqMap.insert(TAADataVT('f', 0.337662337662338 ));
    speerOccurenceFreqMap.insert(TAADataVT('p', 0.480519480519481 ));
    speerOccurenceFreqMap.insert(TAADataVT('s', 0.714285714285714 ));
    speerOccurenceFreqMap.insert(TAADataVT('t', 0.584415584415584 ));
    speerOccurenceFreqMap.insert(TAADataVT('w', 0 ));
    speerOccurenceFreqMap.insert(TAADataVT('y', 0.233766233766234 ));
    speerOccurenceFreqMap.insert(TAADataVT('v', 0.675324675324675 ));
    speerOccurenceFreqMap.insert(TAADataVT('A', 0.818181818181818 ));
    speerOccurenceFreqMap.insert(TAADataVT('R', 0.480519480519481 ));
    speerOccurenceFreqMap.insert(TAADataVT('N', 0.376623376623377 ));
    speerOccurenceFreqMap.insert(TAADataVT('D', 0.493506493506494 ));
    speerOccurenceFreqMap.insert(TAADataVT('C', 0.0779220779220779 ));
    speerOccurenceFreqMap.insert(TAADataVT('Q', 0.350649350649351 ));
    speerOccurenceFreqMap.insert(TAADataVT('E', 0.623376623376623 ));
    speerOccurenceFreqMap.insert(TAADataVT('G', 0.779220779220779 ));
    speerOccurenceFreqMap.insert(TAADataVT('H', 0.116883116883117 ));
    speerOccurenceFreqMap.insert(TAADataVT('I', 0.506493506493506 ));
    speerOccurenceFreqMap.insert(TAADataVT('L', 1 ));
    speerOccurenceFreqMap.insert(TAADataVT('K', 0.584415584415584 ));
    speerOccurenceFreqMap.insert(TAADataVT('M', 0.12987012987013 ));
    speerOccurenceFreqMap.insert(TAADataVT('F', 0.337662337662338 ));
    speerOccurenceFreqMap.insert(TAADataVT('P', 0.480519480519481 ));
    speerOccurenceFreqMap.insert(TAADataVT('S', 0.714285714285714 ));
    speerOccurenceFreqMap.insert(TAADataVT('T', 0.584415584415584 ));
    speerOccurenceFreqMap.insert(TAADataVT('W', 0 ));
    speerOccurenceFreqMap.insert(TAADataVT('Y', 0.233766233766234 ));
    speerOccurenceFreqMap.insert(TAADataVT('V', 0.675324675324675 ));
    speerOccurenceFreqMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerOccurenceFreqMap);

/*
    for (TAADataIt it = speerOccurenceFreqMap.begin(); it != speerOccurenceFreqMap.end(); ++it) {
        cout << "speerOccurenceFreqMap[" << it->first << "] = " << it->second << endl;
    }
*/

}

void SpeerInitializePkA_RCOOH() {

    /*  From ALL-INDEX
pk-a_rcooh
A 0.836555360281195 
R 0.755711775043937 
N 0.639718804920914 
D 1 
C 0.644991212653778 
Q 0.797891036906854 
E 0.963093145869947 
G 0.662565905096661 
H 0.499121265377856 
I 0.845342706502636 
L 0.84182776801406 
K 0.750439367311072 
M 0.746924428822496 
F 0.757469244288225 
P 0 
S 0.67311072056239 
T 0.680140597539543 
W 0.834797891036907 
Y 0.755711775043937 
V 0.854130052724077 
    */

    TAADataMap& speerPkA_RCOOHMap = speerPropertyMap["pk-a_rcooh"];
    speerPkA_RCOOHMap.clear();

    speerPkA_RCOOHMap.insert(TAADataVT('a', 0.836555360281195 ));
    speerPkA_RCOOHMap.insert(TAADataVT('r', 0.755711775043937 ));
    speerPkA_RCOOHMap.insert(TAADataVT('n', 0.639718804920914 ));
    speerPkA_RCOOHMap.insert(TAADataVT('d', 1 ));
    speerPkA_RCOOHMap.insert(TAADataVT('c', 0.644991212653778 ));
    speerPkA_RCOOHMap.insert(TAADataVT('q', 0.797891036906854 ));
    speerPkA_RCOOHMap.insert(TAADataVT('e', 0.963093145869947 ));
    speerPkA_RCOOHMap.insert(TAADataVT('g', 0.662565905096661 ));
    speerPkA_RCOOHMap.insert(TAADataVT('h', 0.499121265377856 ));
    speerPkA_RCOOHMap.insert(TAADataVT('i', 0.845342706502636 ));
    speerPkA_RCOOHMap.insert(TAADataVT('l', 0.84182776801406 ));
    speerPkA_RCOOHMap.insert(TAADataVT('k', 0.750439367311072 ));
    speerPkA_RCOOHMap.insert(TAADataVT('m', 0.746924428822496 ));
    speerPkA_RCOOHMap.insert(TAADataVT('f', 0.757469244288225 ));
    speerPkA_RCOOHMap.insert(TAADataVT('p', 0 ));
    speerPkA_RCOOHMap.insert(TAADataVT('s', 0.67311072056239 ));
    speerPkA_RCOOHMap.insert(TAADataVT('t', 0.680140597539543 ));
    speerPkA_RCOOHMap.insert(TAADataVT('w', 0.834797891036907 ));
    speerPkA_RCOOHMap.insert(TAADataVT('y', 0.755711775043937 ));
    speerPkA_RCOOHMap.insert(TAADataVT('v', 0.854130052724077 ));
    speerPkA_RCOOHMap.insert(TAADataVT('A', 0.836555360281195 ));
    speerPkA_RCOOHMap.insert(TAADataVT('R', 0.755711775043937 ));
    speerPkA_RCOOHMap.insert(TAADataVT('N', 0.639718804920914 ));
    speerPkA_RCOOHMap.insert(TAADataVT('D', 1 ));
    speerPkA_RCOOHMap.insert(TAADataVT('C', 0.644991212653778 ));
    speerPkA_RCOOHMap.insert(TAADataVT('Q', 0.797891036906854 ));
    speerPkA_RCOOHMap.insert(TAADataVT('E', 0.963093145869947 ));
    speerPkA_RCOOHMap.insert(TAADataVT('G', 0.662565905096661 ));
    speerPkA_RCOOHMap.insert(TAADataVT('H', 0.499121265377856 ));
    speerPkA_RCOOHMap.insert(TAADataVT('I', 0.845342706502636 ));
    speerPkA_RCOOHMap.insert(TAADataVT('L', 0.84182776801406 ));
    speerPkA_RCOOHMap.insert(TAADataVT('K', 0.750439367311072 ));
    speerPkA_RCOOHMap.insert(TAADataVT('M', 0.746924428822496 ));
    speerPkA_RCOOHMap.insert(TAADataVT('F', 0.757469244288225 ));
    speerPkA_RCOOHMap.insert(TAADataVT('P', 0 ));
    speerPkA_RCOOHMap.insert(TAADataVT('S', 0.67311072056239 ));
    speerPkA_RCOOHMap.insert(TAADataVT('T', 0.680140597539543 ));
    speerPkA_RCOOHMap.insert(TAADataVT('W', 0.834797891036907 ));
    speerPkA_RCOOHMap.insert(TAADataVT('Y', 0.755711775043937 ));
    speerPkA_RCOOHMap.insert(TAADataVT('V', 0.854130052724077 ));
    speerPkA_RCOOHMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerPkA_RCOOHMap);

/*
    for (TAADataIt it = speerPkA_RCOOHMap.begin(); it != speerPkA_RCOOHMap.end(); ++it) {
        cout << "speerPkA_RCOOHMap[" << it->first << "] = " << it->second << endl;
    }
*/

}
void SpeerInitializePolarity() {

    /*  From ALL-INDEX
polarity
A 0.395061728395062 
R 0.691358024691358 
N 0.827160493827161 
D 1 
C 0.074074074074074 
Q 0.691358024691358 
E 0.91358024691358 
G 0.506172839506173 
H 0.679012345679012 
I 0.037037037037037 
L 0 
K 0.790123456790124 
M 0.0987654320987654 
F 0.037037037037037 
P 0.382716049382716 
S 0.530864197530864 
T 0.45679012345679 
W 0.0617283950617284 
Y 0.160493827160494 
V 0.123456790123457 
    */

    TAADataMap& speerPolarityMap = speerPropertyMap["polarity"];
    speerPolarityMap.clear();

    speerPolarityMap.insert(TAADataVT('a', 0.395061728395062 ));
    speerPolarityMap.insert(TAADataVT('r', 0.691358024691358 ));
    speerPolarityMap.insert(TAADataVT('n', 0.827160493827161 ));
    speerPolarityMap.insert(TAADataVT('d', 1 ));
    speerPolarityMap.insert(TAADataVT('c', 0.074074074074074 ));
    speerPolarityMap.insert(TAADataVT('q', 0.691358024691358 ));
    speerPolarityMap.insert(TAADataVT('e', 0.91358024691358 ));
    speerPolarityMap.insert(TAADataVT('g', 0.506172839506173 ));
    speerPolarityMap.insert(TAADataVT('h', 0.679012345679012 ));
    speerPolarityMap.insert(TAADataVT('i', 0.037037037037037 ));
    speerPolarityMap.insert(TAADataVT('l', 0 ));
    speerPolarityMap.insert(TAADataVT('k', 0.790123456790124 ));
    speerPolarityMap.insert(TAADataVT('m', 0.0987654320987654 ));
    speerPolarityMap.insert(TAADataVT('f', 0.037037037037037 ));
    speerPolarityMap.insert(TAADataVT('p', 0.382716049382716 ));
    speerPolarityMap.insert(TAADataVT('s', 0.530864197530864 ));
    speerPolarityMap.insert(TAADataVT('t', 0.45679012345679 ));
    speerPolarityMap.insert(TAADataVT('w', 0.0617283950617284 ));
    speerPolarityMap.insert(TAADataVT('y', 0.160493827160494 ));
    speerPolarityMap.insert(TAADataVT('v', 0.123456790123457 ));
    speerPolarityMap.insert(TAADataVT('A', 0.395061728395062 ));
    speerPolarityMap.insert(TAADataVT('R', 0.691358024691358 ));
    speerPolarityMap.insert(TAADataVT('N', 0.827160493827161 ));
    speerPolarityMap.insert(TAADataVT('D', 1 ));
    speerPolarityMap.insert(TAADataVT('C', 0.074074074074074 ));
    speerPolarityMap.insert(TAADataVT('Q', 0.691358024691358 ));
    speerPolarityMap.insert(TAADataVT('E', 0.91358024691358 ));
    speerPolarityMap.insert(TAADataVT('G', 0.506172839506173 ));
    speerPolarityMap.insert(TAADataVT('H', 0.679012345679012 ));
    speerPolarityMap.insert(TAADataVT('I', 0.037037037037037 ));
    speerPolarityMap.insert(TAADataVT('L', 0 ));
    speerPolarityMap.insert(TAADataVT('K', 0.790123456790124 ));
    speerPolarityMap.insert(TAADataVT('M', 0.0987654320987654 ));
    speerPolarityMap.insert(TAADataVT('F', 0.037037037037037 ));
    speerPolarityMap.insert(TAADataVT('P', 0.382716049382716 ));
    speerPolarityMap.insert(TAADataVT('S', 0.530864197530864 ));
    speerPolarityMap.insert(TAADataVT('T', 0.45679012345679 ));
    speerPolarityMap.insert(TAADataVT('W', 0.0617283950617284 ));
    speerPolarityMap.insert(TAADataVT('Y', 0.160493827160494 ));
    speerPolarityMap.insert(TAADataVT('V', 0.123456790123457 ));
    speerPolarityMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerPolarityMap);

/*
    for (TAADataIt it = speerPolarityMap.begin(); it != speerPolarityMap.end(); ++it) {
        cout << "speerPolarityMap[" << it->first << "] = " << it->second << endl;
    }
*/

}
void SpeerInitializeSize() {

    /*  From ALL-INDEX
size
A 0.285714285714286 
R 1 
N 0.642857142857143 
D 0.285714285714286 
C 0.357142857142857 
Q 0.785714285714286 
E 0.642857142857143 
G 0 
H 0.785714285714286 
I 0.714285714285714 
L 0.714285714285714 
K 0.928571428571428 
M 0.785714285714286 
F 0.857142857142857 
P 0.714285714285714 
S 0.357142857142857 
T 0.642857142857143 
W 0.928571428571428 
Y 0.928571428571428 
V 0.642857142857143 
    */

    TAADataMap& speerSizeMap = speerPropertyMap["size"];
    speerSizeMap.clear();

    speerSizeMap.insert(TAADataVT('a', 0.285714285714286 ));
    speerSizeMap.insert(TAADataVT('r', 1 ));
    speerSizeMap.insert(TAADataVT('n', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('d', 0.285714285714286 ));
    speerSizeMap.insert(TAADataVT('c', 0.357142857142857 ));
    speerSizeMap.insert(TAADataVT('q', 0.785714285714286 ));
    speerSizeMap.insert(TAADataVT('e', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('g', 0 ));
    speerSizeMap.insert(TAADataVT('h', 0.785714285714286 ));
    speerSizeMap.insert(TAADataVT('i', 0.714285714285714 ));
    speerSizeMap.insert(TAADataVT('l', 0.714285714285714 ));
    speerSizeMap.insert(TAADataVT('k', 0.928571428571428 ));
    speerSizeMap.insert(TAADataVT('m', 0.785714285714286 ));
    speerSizeMap.insert(TAADataVT('f', 0.857142857142857 ));
    speerSizeMap.insert(TAADataVT('p', 0.714285714285714 ));
    speerSizeMap.insert(TAADataVT('s', 0.357142857142857 ));
    speerSizeMap.insert(TAADataVT('t', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('w', 0.928571428571428 ));
    speerSizeMap.insert(TAADataVT('y', 0.928571428571428 ));
    speerSizeMap.insert(TAADataVT('v', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('A', 0.285714285714286 ));
    speerSizeMap.insert(TAADataVT('R', 1 ));
    speerSizeMap.insert(TAADataVT('N', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('D', 0.285714285714286 ));
    speerSizeMap.insert(TAADataVT('C', 0.357142857142857 ));
    speerSizeMap.insert(TAADataVT('Q', 0.785714285714286 ));
    speerSizeMap.insert(TAADataVT('E', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('G', 0 ));
    speerSizeMap.insert(TAADataVT('H', 0.785714285714286 ));
    speerSizeMap.insert(TAADataVT('I', 0.714285714285714 ));
    speerSizeMap.insert(TAADataVT('L', 0.714285714285714 ));
    speerSizeMap.insert(TAADataVT('K', 0.928571428571428 ));
    speerSizeMap.insert(TAADataVT('M', 0.785714285714286 ));
    speerSizeMap.insert(TAADataVT('F', 0.857142857142857 ));
    speerSizeMap.insert(TAADataVT('P', 0.714285714285714 ));
    speerSizeMap.insert(TAADataVT('S', 0.357142857142857 ));
    speerSizeMap.insert(TAADataVT('T', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('W', 0.928571428571428 ));
    speerSizeMap.insert(TAADataVT('Y', 0.928571428571428 ));
    speerSizeMap.insert(TAADataVT('V', 0.642857142857143 ));
    speerSizeMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerSizeMap);

/*
    for (TAADataIt it = speerSizeMap.begin(); it != speerSizeMap.end(); ++it) {
        cout << "speerSizeMap[" << it->first << "] = " << it->second << endl;
    }
*/

}

void SpeerInitializeVolume() {

    /*  From ALL-INDEX
volume
A 0.167664670658683 
R 0.724550898203593 
N 0.317365269461078 
D 0.305389221556886 
C 0.311377245508982 
Q 0.491017964071856 
E 0.479041916167665 
G 0 
H 0.55688622754491 
I 0.646706586826347 
L 0.646706586826347 
K 0.694610778443114 
M 0.610778443113772 
F 0.772455089820359 
P 0.176646706586826 
S 0.173652694610778 
T 0.347305389221557 
W 1 
Y 0.796407185628743 
V 0.485029940119761 
    */

    TAADataMap& speerVolumeMap = speerPropertyMap["volume"];
    speerVolumeMap.clear();

    speerVolumeMap.insert(TAADataVT('a', 0.167664670658683 ));
    speerVolumeMap.insert(TAADataVT('r', 0.724550898203593 ));
    speerVolumeMap.insert(TAADataVT('n', 0.317365269461078 ));
    speerVolumeMap.insert(TAADataVT('d', 0.305389221556886 ));
    speerVolumeMap.insert(TAADataVT('c', 0.311377245508982 ));
    speerVolumeMap.insert(TAADataVT('q', 0.491017964071856 ));
    speerVolumeMap.insert(TAADataVT('e', 0.479041916167665 ));
    speerVolumeMap.insert(TAADataVT('g', 0 ));
    speerVolumeMap.insert(TAADataVT('h', 0.55688622754491 ));
    speerVolumeMap.insert(TAADataVT('i', 0.646706586826347 ));
    speerVolumeMap.insert(TAADataVT('l', 0.646706586826347 ));
    speerVolumeMap.insert(TAADataVT('k', 0.694610778443114 ));
    speerVolumeMap.insert(TAADataVT('m', 0.610778443113772 ));
    speerVolumeMap.insert(TAADataVT('f', 0.772455089820359 ));
    speerVolumeMap.insert(TAADataVT('p', 0.176646706586826 ));
    speerVolumeMap.insert(TAADataVT('s', 0.173652694610778 ));
    speerVolumeMap.insert(TAADataVT('t', 0.347305389221557 ));
    speerVolumeMap.insert(TAADataVT('w', 1 ));
    speerVolumeMap.insert(TAADataVT('y', 0.796407185628743 ));
    speerVolumeMap.insert(TAADataVT('v', 0.485029940119761 ));
    speerVolumeMap.insert(TAADataVT('A', 0.167664670658683 ));
    speerVolumeMap.insert(TAADataVT('R', 0.724550898203593 ));
    speerVolumeMap.insert(TAADataVT('N', 0.317365269461078 ));
    speerVolumeMap.insert(TAADataVT('D', 0.305389221556886 ));
    speerVolumeMap.insert(TAADataVT('C', 0.311377245508982 ));
    speerVolumeMap.insert(TAADataVT('Q', 0.491017964071856 ));
    speerVolumeMap.insert(TAADataVT('E', 0.479041916167665 ));
    speerVolumeMap.insert(TAADataVT('G', 0 ));
    speerVolumeMap.insert(TAADataVT('H', 0.55688622754491 ));
    speerVolumeMap.insert(TAADataVT('I', 0.646706586826347 ));
    speerVolumeMap.insert(TAADataVT('L', 0.646706586826347 ));
    speerVolumeMap.insert(TAADataVT('K', 0.694610778443114 ));
    speerVolumeMap.insert(TAADataVT('M', 0.610778443113772 ));
    speerVolumeMap.insert(TAADataVT('F', 0.772455089820359 ));
    speerVolumeMap.insert(TAADataVT('P', 0.176646706586826 ));
    speerVolumeMap.insert(TAADataVT('S', 0.173652694610778 ));
    speerVolumeMap.insert(TAADataVT('T', 0.347305389221557 ));
    speerVolumeMap.insert(TAADataVT('W', 1 ));
    speerVolumeMap.insert(TAADataVT('Y', 0.796407185628743 ));
    speerVolumeMap.insert(TAADataVT('V', 0.485029940119761 ));
    speerVolumeMap.insert(TAADataVT('-', -1));  //  gap


    SpeerComputeBZXPropertyValues(speerVolumeMap);

/*
    for (TAADataIt it = speerVolumeMap.begin(); it != speerVolumeMap.end(); ++it) {
        cout << "speerVolumeMap[" << it->first << "] = " << it->second << endl;
    }
*/

}


//  Data needed for evolutionary-distance scoring
void SpeerInitializeDefaultEDistData() {
    SpeerInitializeContactNumber();
    SpeerInitializeAbsoluteEntropy();
    SpeerInitializeAccSurfArea();
    SpeerInitializeHydrationPotl();
    SpeerInitializeHydropathyIndex();
    SpeerInitializeHydrophobicityIndex();
    SpeerInitializeIsoElecPoint();
    SpeerInitializeMolWeight();
    SpeerInitializeNetCharge();
    SpeerInitializeNormFlexParam();
    SpeerInitializeRelativeMutability();
    SpeerInitializeAccSurfAreaInFoldedProtein();
    SpeerInitializeSideChainOrientPref();
    SpeerInitializeOccurenceFreq();
    SpeerInitializePkA_RCOOH();
    SpeerInitializePolarity();
    SpeerInitializeSize();
    SpeerInitializeVolume();
}


void SpeerInitializePairDistanceMap()
{
    unsigned int i, j;
    TSpeerScore sum, term;
    TPropertyMapIt propIt, propEnd = speerPropertyMap.end();
    string::value_type iChar, jChar;

    for (i = 0; i < nAA; ++i) {
        iChar = AminoAcidChar[i];
        speerPairDistanceMap.insert(TPairDistanceMapVT(TResiduePair(iChar, iChar), 0.0));

        for (j = 0; j < i; ++j) {
            sum = 0.0;
            propIt = speerPropertyMap.begin();
            jChar = AminoAcidChar[j];

            for (; propIt != propEnd; ++propIt) {
                TAADataMap& dataMap = propIt->second;
                term = dataMap[iChar] - dataMap[jChar];
                sum += term*term;
            }
            sum = sqrt(sum);
            speerPairDistanceMap.insert(TPairDistanceMapVT(TResiduePair(iChar, jChar), sum));
            speerPairDistanceMap.insert(TPairDistanceMapVT(TResiduePair(jChar, iChar), sum));

//            cerr << iChar << ", " << jChar << ":  " << sum<< endl;
        }
    }
}


void SpeerInitializeDefaultGlobalData() {

    SpeerInitializeBackgroundFreqData();
    SpeerInitializeDefaultEDistData();

    SpeerInitializePairDistanceMap();
}


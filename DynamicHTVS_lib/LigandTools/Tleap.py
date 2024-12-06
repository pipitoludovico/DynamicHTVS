from subprocess import run, CalledProcessError, DEVNULL


def TleapLigand(mol2path: str, name: str) -> None:
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",

         "loadamberparams UNL.frcmod",
         f"UNL = loadmol2 {mol2path}",
         "check UNL",
         "saveoff UNL UNL.lib",
         "saveamberparm UNL ./gbsa/ligand.prmtop ./gbsa/ligand.inpcrd",
         "savepdb UNL ./gbsa/ligand.pdb",
         "quit"]
    WriteTleap(_, name)


def TleapReceptor(recPdbPath: str, name: str) -> None:
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",

         f'rec = loadpdb {recPdbPath}',
         "saveamberparm rec ./gbsa/receptor.prmtop ./gbsa/receptor.inpcrd",
         "savepdb rec ./gbsa/ligand.pdb", "quit"]
    WriteTleap(_, name)


def TleapMakeGBSAComplex(recPdbPath: str, mol2path: str, name: str) -> None:
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "loadamberparams UNL.frcmod",
         "set default PBRadii mbondi3",
         f"UNL = loadmol2 {mol2path}",
         "check UNL",
         f"rec = loadpdb {recPdbPath}",
         "complex = combine{rec UNL}",
         "saveamberparm complex ./gbsa/complex.prmtop ./gbsa/complex.inpcrd",
         "savepdb complex ./gbsa/complex.pdb",
         "quit"]
    WriteTleap(_, name)


def TleapMakeComplexCLASH(recPdbPath: str, mol2path: str, name: str) -> None:
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",
         "loadamberparams UNL.frcmod",
         f"UNL = loadmol2 {mol2path}",
         "check UNL",
         f"rec = loadpdb {recPdbPath}",
         "complex = combine{rec UNL}",
         'setBox complex "vdw"',
         'setBox complex "vdw"',
         "savepdb complex ./clashed.pdb", "quit"]
    WriteTleap(_, name)


def TleapMakeComplexMD(complexNOCLASHpath: str, mol2path: str, name: str) -> None:
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",
         "loadamberparams UNL.frcmod",
         f"UNL = loadmol2 {mol2path}",
         "check UNL",
         f"complex = loadpdb {complexNOCLASHpath}",
         'setBox complex "vdw"',
         "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd",
         "quit"]
    WriteTleap(_, name)


def TleapIonize(mol2path, complexPdbPath, conc, name):
    """Deprecated"""
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",
         "loadamberparams UNL.frcmod",
         f"UNL = loadmol2 {mol2path}",
         "check UNL",
         f"complex = loadpdb {complexPdbPath}",
         "check complex",
         "saveamberparm complex ./gbsa/complex.prmtop ./gbsa/complex.inpcrd",
         # now we save for cMD
         'setBox complex "vdw"',
         'setBox complex "vdw"',
         "solvateBox complex TIP3PBOX 15 0.75",
         "addIons complex Na+ 0",
         "addIons complex Cl- 0",
         f"addIons complex Na+ {conc}",
         f"addIons complex Cl- {conc}",
         "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd",
         "savepdb complex ./system/complex.pdb"]
    WriteTleap(_, name)


def WriteTleap(_: list, name):
    with open(f'inleap_{name}', 'w') as inleap:
        for tleapCommand in _:
            inleap.write(tleapCommand + "\n")
    try:
        run(f"tleap -f inleap_{name}", shell=True, stdout=DEVNULL)
    except CalledProcessError:
        print("tleap failed")
        exit()

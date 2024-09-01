from os import makedirs
from subprocess import run, CalledProcessError, DEVNULL


def RunTleap(structureFile, solvate, ionize=None, conc=1, MOL2=None) -> None:
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",

         "loadamberparams UNL.frcmod" if structureFile.endswith('.mol2') or "complex" in structureFile else "",
         f"UNL = loadmol2 {MOL2}" if MOL2 else "",
         f"UNL = loadmol2 {structureFile}" if structureFile.endswith('.mol2') and "complex" not in structureFile else "",
         "check UNL" if structureFile.endswith('.mol2') and "complex" not in structureFile else "",
         "saveoff UNL UNL.lib" if structureFile.endswith('.mol2') and "complex" not in structureFile else "",
         "saveamberparm UNL ./gbsa/ligand.prmtop ./gbsa/ligand.inpcrd" if structureFile.endswith('.mol2') and "complex" not in structureFile else "",
         "savepdb UNL ./gbsa/ligand.pdb" if structureFile.endswith('.mol2') and "complex" not in structureFile else "",

         'rec = loadpdb ../../receptor/system.pdb' if structureFile.endswith('.pdb') and "system" in structureFile else "",
         "saveamberparm rec ./gbsa/receptor.prmtop ./gbsa/receptor.inpcrd" if structureFile.endswith('.pdb') and "system" in structureFile else "",
         "savepdb rec ./gbsa/receptor.pdb" if structureFile.endswith('.pdb') and "system" in structureFile else "",

         "complex = loadpdb complex.pdb" if structureFile.endswith('.pdb') and "complex" in structureFile else "",
         "saveamberparm complex ./gbsa/complex.prmtop ./gbsa/complex.inpcrd" if structureFile.endswith('.pdb') and not ionize and "complex" in structureFile else "",
         'setBox complex "vdw"' if "complex" in structureFile else "",
         'setBox complex "vdw"' if solvate else "",
         "solvateBox complex TIP3PBOX 15 0.75" if solvate else "",
         "savepdb complex solvated.pdb" if ionize and "complex" in structureFile else "",

         "complex = loadpdb solvated.pdb" if ionize and "complex" in structureFile else "",
         "addIons complex Na+ 0" if ionize and "complex" in structureFile else "",
         "addIons complex Cl- 0" if ionize and "complex" in structureFile else "",
         f"addIons complex Na+ {conc}" if ionize and "complex" in structureFile else "",
         f"addIons complex Cl- {conc}" if ionize and "complex" in structureFile else "",
         'setBox complex "vdw"' if ionize and "complex" in structureFile else "",
         "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd" if ionize and "complex" in structureFile else "",
         "savepdb complex ./system/complex.pdb" if ionize and "complex" in structureFile else "",
         "quit"]

    with open('inleap', 'w') as inleap:
        for tleapCommand in _:
            inleap.write(tleapCommand + "\n")
    try:
        run("tleap -f inleap", shell=True, stdout=DEVNULL)
    except CalledProcessError:
        print("tleap failed")
        exit()

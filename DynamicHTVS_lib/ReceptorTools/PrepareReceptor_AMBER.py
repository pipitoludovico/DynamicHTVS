import os.path

from subprocess import run, DEVNULL


def PrepareProtein() -> None:
    """Prepares the system for AMBER using htmd."""
    if not os.path.exists('./receptor/.check'):
        run("sed -i 's/HSP/HIS/g' ./receptor/system.pdb", shell=True)
        run("sed -i 's/HSD/HIS/g' ./receptor/system.pdb", shell=True)
        run("sed -i 's/HSE/HIS/g' ./receptor/system.pdb", shell=True)
        import htmd.ui as htmdmodule
        print("Building the receptor from MD...\n\n")
        tick = None
        with open("./receptor/system.pdb", 'r') as originalREC:
            for line in originalREC:
                if 'MEMB' in line:
                    tick = 32.0
        try:
            protein = htmdmodule.Molecule('./receptor/system.pdb')
            receptor = htmdmodule.systemPrepare(protein, hydrophobic_thickness=tick, ignore_ns_errors=True,
                                                hold_nonpeptidic_bonds=True, titration=True)
            receptor.write('./receptor/system_H.pdb')
            # we need to remove hydrogens as systemPrepare adds CHARMM-like Hs atomtypes to the system...
            run('pdb4amber -i ./receptor/system_H.pdb -o ./receptor/system.pdb -y; touch ./receptor/.check',
                shell=True, stdout=DEVNULL)
            del protein, receptor  # clearing memory as we don't need those anymore
        except Exception as e:
            print("Building the system raised this exception:\n\n", e)
            print("\nTrying pdb4amber first.\n")
            run('pdb4amber -i ./receptor/system.pdb -o ./receptor/receptor_amber.pdb -y -a', shell=True, stdout=DEVNULL)
            protein = htmdmodule.Molecule('./receptor/receptor_amber.pdb')
            receptor = htmdmodule.systemPrepare(protein, hydrophobic_thickness=tick, ignore_ns_errors=True, hold_nonpeptidic_bonds=True, titration=True)
            receptor.write('./receptor/receptor_H.pdb')
            run('pdb4amber -i ./receptor/receptor_H.pdb -o ./receptor/system.pdb -y; touch ./receptor/.check',
                shell=True)
    else:
        print("Found hidden .check file from a previous pdb4amber run. Remove this file if you want to run pdb4amber again")

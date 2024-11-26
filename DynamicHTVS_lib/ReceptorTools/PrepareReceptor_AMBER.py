import os
from os import path, listdir, makedirs

from subprocess import run, DEVNULL


def PrepareProtein() -> None:
    """Prepares the system for AMBER using pdb2pqr."""
    makedirs('./receptor/backups', exist_ok=True)
    pdbFilelist: list = []
    if not path.exists('./receptor/.check'):
        if any(file.endswith('psf') for file in os.listdir('./receptor')):
            print("Found a CHARMM psf files while -amber keyword is on. Attempting to convert them to AMBER.")

            for file in listdir('./receptor'):
                if file.endswith('pdb'):
                    run(f"cp ./receptor/{file} ./receptor/backups/{file}", shell=True)
                    run(f"sed -i 's/HSP/HIS/g' ./receptor/{file}", shell=True)
                    run(f"sed -i 's/HSD/HIS/g' ./receptor/{file}", shell=True)
                    run(f"sed -i 's/HSE/HIS/g' ./receptor/{file}", shell=True)
                    run(f"sed -i 's/HIP/HIS/g' ./receptor/{file}", shell=True)
                    pdbFilelist.append(f"./receptor/{file}")
        for file in pdbFilelist:
            try:
                with open('./logs/pdb4amberlog.log', 'a') as pdb4amberlogger:
                    run(f'pdb4amber -i {file} -o ./receptor/tmp.pdb -a;', shell=True, stdout=DEVNULL, stderr=pdb4amberlogger)
                    run('pdb2pqr ./receptor/tmp.pdb ./receptor/tmp.pqr --ff AMBER --ffout AMBER --titration-state-method propka', shell=True, stdout=DEVNULL, stderr=pdb4amberlogger)
                    run(f'cpptraj  -p ./receptor/tmp.pqr -y ./receptor/tmp.pqr -x ./receptor/tmp.pdb;mv ./receptor/tmp.pdb {file}', shell=True, stdout=DEVNULL, stderr=pdb4amberlogger)
                    print(f"\nConversion of {file} complete. Check your system before proceding.")
            except Exception as e:
                print(e)
        run("touch ./receptor/.check", shell=True)
    else:
        print("Found hidden .check file from a previous pdb4amber run. Remove this file if you want to run pdb4amber again")

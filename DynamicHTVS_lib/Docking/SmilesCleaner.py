from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover


class SmilesCleaner:
    def __init__(self, smiles):
        self.smiles = smiles
        self.molObjs = []
        self.canonicalSmiles = []
        self.remover = SaltRemover()

    def getCanonicalSmiles(self):
        try:
            self.molObjs = [Chem.MolFromSmiles(smi) for smi in self.smiles]
            for mol in self.molObjs:
                Chem.RemoveStereochemistry(mol)
                self.remover.StripMol(mol)
            self.canonicalSmiles = [Chem.MolToSmiles(mol) for mol in self.molObjs]

        except Exception as e:
            print(e)
            print("\n\n********** Check your SMILES or pdb file **********\n\n")
            exit()
        return self.canonicalSmiles

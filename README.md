# DockTBuilder

split_molecules -- functions

- pdb_exist(pathName): verifica se o arquivo PDB existe no caminho especificado
   
- split_file(pathName): divide o arquivo PDB em dois dataFrames, um para átomos e outro para hetátomos agrupados por nome do resíduo
    
- add_hydrogens(hetatms_df): adiciona hidrogênios a um dataFrame de hetero átomos e retorna o dataFrame resultante
    
- save_mol_as_pdb(mol, output_file): salva uma molécula (objeto Mol) em um arquivo PDB especificado
    
- add_hydrogens(hetatms_df): adiciona hidrogênios a um dataFrame de hetátomos, retornando o Mol com hidrogênios 
  

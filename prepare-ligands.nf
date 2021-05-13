params.ligands = 'ligands.smi'
params.chunk_size = 5000
params.min_hac = 0
params.max_hac = 100000
params.min_ph = 5.0
params.max_ph = 9.0
params.num_charges = 1

smilesfile = file(params.ligands)


process splitter {

    container 'informaticsmatters/obabel:latest'
    publishDir './ligands', mode: 'copy'

    input:
    file smiles from smilesfile

    output:
    file 'x*.smi' into chunks

    """
    split -l $params.chunk_size -d -a 6 --additional-suffix .smi $smiles
    """
}

process enumerate {

    container 'informaticsmatters/rdkit_pipelines:latest'
    publishDir './ligands', mode: 'copy'

    input:
    file chunks from chunks.flatten()

    output:
    file 'enumerated_*.smi' optional true into enumerated

    """
    FILE=enumerated_${chunks.name[0..-5]}.smi
    python -m pipelines.rdkit.enumerate_candidates -i '$chunks' -o \$FILE\
      --enumerate-charges --enumerate-chirals --enumerate-tautomers --name-column 1 --num-charges $params.num_charges\
      --min-hac ${params.min_hac} --max-hac ${params.max_hac}\
      --min-ph ${params.min_ph} --max-ph ${params.max_ph}
    if [ ! -s \$FILE ]
    then
      echo "File is empty - deleting"
      rm \$FILE
    fi  
    """
}


process prepareObabel {
    
    container 'informaticsmatters/obabel:latest'
    publishDir './ligands', mode: 'move'

    input:
    file molecules from enumerated.flatten()

    output:
    file 'Prep_*.sdf' into prepared_candidates

    """
    obabel '$molecules' -Omols.sdf --gen2D > obabel1.log
    obabel mols.sdf -O'Prep_${molecules.name[0..-5]}.sdf' --gen3D > obabel2.log
    """
}


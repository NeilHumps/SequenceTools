import os
import shutil
import subprocess
import tempfile

from seqtools.general import load_csv

GENE_FILE = os.path.join(
    os.split(__file__)[0],
    '../data/gene.txt',
)

WORK_DIR = tempfile.mkdtemp()

QSUB_HEADERS = '#!/bin/bash\n' \
               '#$ -N "trinity"\n' \
               '#$ -e errors.log\n' \
               '#$ -cwd\n\n' \
               'module load trinity/2.1.1\n' \
               'module load samtools/1.3\n\n'

TRINITY_CMD = 'Trinity ' \
              '--genome_guided_bam {0} ' \
              '--genome_guided_max_intron 100000 ' \
              '--CPU 8 --max_memory 20G ' \
              '--output {1}\n'


def run_trinity_gene_qsub(gene_name, bam_file, output_dir,
                          pad=1000, # expand gene coordinates by 1kb
                          gene_file=GENE_FILE, work_dir=WORK_DIR,
                          qsub_headers=QSUB_HEADERS, trinity_cmd=TRINITY_CMD):
    _, genes = load_csv(gene_file, sep='\t')
    gene = [i for i in genes if i['gene_name'] == gene_name][0]
    gene_pos = '{0}:{1}-{2}'.format(
      gene['chr'], int(gene['start'])-int(pad), int(gene['end'])+int(pad)
    )
    qsub_file = os.path.join(work_dir, 'trinity.sh')
    bam_tmp = os.path.join(work_dir, gene_name + '.bam')
    with open(qsub_file, 'w') as output_fh:
      output_fh.write(QSUB_HEADERS)
      output_fh.write('samtools view -h -b {0} {1} > {2}\n'.format(
          bam_file, gene_pos, bam_tmp,
      ))
      output_fh.write(TRINITY_CMD.format(bam_tmp, output_dir))
      #subprocess.check_call('echo hello', shell=True)
      subprocess.check_call(['qsub', qsub_file], shell=True)
      # except subprocess.CalledProcessError:
      #     import pdb; pdb.set_trace()
      shutil.rmtree(work_dir)


if __name__ == '__main__':
    pass

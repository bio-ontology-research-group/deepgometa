import click as ck

@ck.command()
def main():
    ncbi_divs = set(['0', '3', '9'])
    with open('data/taxdump/nodes.dmp') as f:
        for line in f:
            it = [x.strip() for x in line.split('|')]
            if it[4] in ncbi_divs:
                print(it[0])



if __name__ == '__main__':
    main()

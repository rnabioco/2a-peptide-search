#!/usr/bin/env python3
"""
Merge multiple Stockholm alignments into a single alignment.
Uses esl-alimerge under the hood but provides a cleaner interface.
"""

import subprocess
import tempfile
import click


@click.command()
@click.argument('alignments', nargs=-1, type=click.Path(exists=True))
@click.argument('output', type=click.Path())
def main(alignments, output):
    """Merge multiple alignments into one."""

    # Create temporary file with alignment list
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        for aln in alignments:
            tmp.write(f"{aln}\n")
        tmp_name = tmp.name

    # Run esl-alimerge
    try:
        with open(output, 'w') as out:
            subprocess.run(
                ['esl-alimerge', '--list', tmp_name],
                stdout=out,
                check=True
            )
        click.echo(f"Merged {len(alignments)} alignments -> {output}")
    finally:
        import os
        os.unlink(tmp_name)


if __name__ == '__main__':
    main()

def convert_fasta(handle):
    """
    Import FASTA-formatted sequences from file stream
    :param handle:  File stream in read mode, or the contents of a file split into lines
    :return: list of lists containing header-sequence pairs
    """
    result = []
    sequence, h = '', ''

    # Verifies files have correct formatting
    for i, line in enumerate(handle):
        if line.startswith('$'):
            continue
        elif line.startswith('>') or line.startswith('#'):
            break
        else:
            print("No header")
            raise NameError

    if hasattr(handle, 'seek'):
        handle.seek(0)

    for line in handle:
        if line.startswith('$'):  # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''  # reset
            h = line.strip('>#\n\r')
        else:
            sequence += line.strip('\n\r').upper()

    result.append([h, sequence])  # handle last entry
    return result

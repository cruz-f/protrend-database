def args_length(*args):
    size = 0

    for param in args:

        if param is not None:
            param_size = len(param)
            if param_size > size:
                size = param_size

    return size


def scale_arg(arg, size):
    if arg is None:
        return [None] * size

    elif len(arg) == 0:
        return [None] * size

    elif len(arg) == 1:
        return arg * size

    elif len(arg) != size:
        raise ValueError(f'Invalid input size of {len(arg)}')

    return arg

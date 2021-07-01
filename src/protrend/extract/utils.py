def parsing_spider_arguments(argument: str):

    if not argument:
        return

    children = argument.split(',')

    return tuple(children)

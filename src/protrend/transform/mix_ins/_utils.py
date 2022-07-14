def get_values(df, col):
    series = df.get(col)

    if series is None:
        return

    return series.to_list()

def allSame(strand):
    assert len(strand) > 0
    n = strand[0]
    for i in strand:
        if i != n:
            return False
    return True
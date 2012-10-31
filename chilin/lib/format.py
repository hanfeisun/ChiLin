def group (source, n):
    """
    group([1,2,3,4],2) -> [[1, 2], [3, 4]]
    group([1,2,3,4],3) -> [[1, 2, 3], [4]]
    """
    acc = []
    rest = source
    if not n:
        return []
    while rest:
        acc.append(source[:n])
        rest = source[n:]
        source = rest        
    return acc

"""
digit_format(1234567890) -> '1,234,567,890'
digit_format(123456789) -> '123,456,789'
digit_format(12) -> '12'
digit_format(21124324,"\,") -> '21\\,124\\,324'
"""
digit_format = lambda num,sep=",": sep[::-1].join(group(str(num)[::-1],3))[::-1]



"""percent_format(0.8964) -> '89.64%'"""
percent_format = lambda num: str(num*100)+"%"

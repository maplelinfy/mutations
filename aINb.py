
#arr should be the form [[chr, pos], [chr, pos]...]
def a_in_b(arr_a, arr_b):
    arr_a = sorted(arr_a, key=lambda x: (x[0], x[1]))
    arr_b = sorted(arr_b, key=lambda x: (x[0], x[1]))
    count = 0
    index = 0
    for i in range(len(arr_a)):
        chrr = arr_a[i][0]
        pos = arr_a[i][1]
        while (chrr == arr_b[index][0] and pos > arr_b[index][1]) or (chrr > arr_b[index][0]):
            index += 1
            if index == len(arr_b):
                break
        if index < len(arr_b):
            if chrr != arr_b[index][0] or pos != arr_b[index][1]:
                continue
            else:
                count += 1
        else:
            break
    print(count)

# arr_a should be the form[[chr, pos], [chr, pos]...]
# arr_b should be the form[[chr, start, end], [chr, start, end]...]
def a_in_b_region(arr_a, arr_b):
    arr_a = sorted(arr_a, key=lambda x: (x[0], x[1]))
    arr_b = sorted(arr_b, key=lambda x: (x[0], x[1]))
    index = 0
    arr = []
    for i in range(len(arr_a)):
        chrr = arr_a[i][0]
        pos = arr_a[i][1]
        while (chrr == arr_b[index][0] and pos > arr_b[index][2]) or (chrr > arr_b[index][0]):
            index += 1
            if index == len(arr_b):
                break
        if index < len(arr_b):
            if chrr != arr_b[index][0] or pos < arr_b[index][1]:
                continue
            else:
                arr.append(arr_a[i])
        else:
            break
    return arr

#arr should be the form [[chr, start, end], [chr, start, end]...]
def calUnionArea(arr_a, arr_b):
    area = 0
    index = 0
    for i in range(len(arr_a)):
        while (arr_b[index][0] < arr_a[i][0]) or (arr_b[index][2] < arr_a[i][1] and arr_b[index][0] == arr_a[i][0]):
            index += 1
            if index == len(arr_b):
                break
        if index < len(arr_b):
            if (arr_b[index][2] - arr_a[i][1]) * (arr_a[i][2] - arr_b[index][1]) >= 0 and arr_b[index][0] == arr_a[i][0]:
                area += min(arr_a[i][2], arr_b[index][2]) - max(arr_a[i][1], arr_b[index][1])
        else:
            break
    return area
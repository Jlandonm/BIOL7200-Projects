
def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:

    horizontal = seq_a
    vertical = seq_b

    #Running helper to fill matrix
    matrix = fill_matrix(horizontal, vertical, match, mismatch, gap)

    #Running helper function to backtrace
    h_align, v_align, score = backtrace(matrix, horizontal, vertical)

    return ((h_align, v_align), score)

def fill_matrix(horizontal, vertical, match, mismatch, gap): 
    rows, cols = (len(vertical) + 1, len(horizontal) + 1)
    matrix = [[[0] for i in range(cols)] for j in range(rows)]

    #Initialize, doing the gaps first
    for i in range(1,rows):
        matrix[i][0][0] = matrix[i - 1][0][0] + gap
    
    for j in range(1,cols):
        matrix[0][j][0] = matrix[0][j - 1][0] + gap

    #Now start the algorithm

    for j in range(1, rows):
        for i in range(1, cols):
            dict = {}
            square = matrix[j][i]
            dict['left'] = matrix[j][i - 1][0] + gap
            dict['up'] = matrix[j - 1][i][0] + gap
            #match
            if horizontal[i - 1] == vertical[j - 1]:
                dict['diagonal'] = matrix[j - 1][i - 1][0] + match
            else:
                dict['diagonal'] = matrix[j - 1][i - 1][0] + mismatch
            #Get order of movement and the score values
            order = [key for (key, value) in sorted(dict.items(), key=lambda key_value: key_value[1], reverse=True)]
            values = [value for (key, value) in sorted(dict.items(), key=lambda key_value: key_value[1], reverse=True)]

            #Update matrix value
            square[0] = values[0]

            #Add this to our square
            square.append(order)
            square.append(values)
    
    return matrix

def backtrace(matrix, horizontal, vertical):

    h_align = ""
    v_align = ""
    ind_h = len(horizontal)
    ind_v = len(vertical)
    score = matrix[len(vertical)][len(horizontal)][0]
    #Now backtrace to get alignment
    while (ind_h > 0 and ind_v > 0):
        #GO
        square = matrix[ind_v][ind_h]
        #Now we need to go through the arrows. Our default is going to be a diagonal > left > up
        #If we only have one arrow -->
        if square[2][0] > square[2][1]:
            move = square[1][0]
        #Two way tie
        elif square[2][0] == square[2][1] and square[2][0] > square[2][2]:
            if square[1][2] == "diagonal":
                move = "left"
            else:
                move = "diagonal"
        #First and third equal, three way tie
        else:
            move = "diagonal"
        #Now do the move and record the base
        if move == "diagonal":
            h_align += horizontal[ind_h - 1]
            v_align += vertical[ind_v - 1]
            ind_h -= 1
            ind_v -= 1
        elif move == "left":
            h_align += horizontal[ind_h - 1]
            v_align += "-"
            ind_h -= 1
        else:
            h_align += "-"
            v_align += vertical[ind_v - 1]
            ind_v -= 1

    #Now we need to reverse the sequences
    h_align = reversed(h_align)
    h_align = "".join(h_align)
    v_align = reversed(v_align)
    v_align = "".join(v_align)


    return h_align, v_align, score
def find_repeats(string, min_length=4, max_length=20):
    '''
    Function to find repeated substrings in a string and return them as a dictionary with number of repeats.
    The substrings must be between min_length and max_length characters long.
    '''

    substrings = [string[i:i+j] for j in range(min_length, max_length) for i in range(len(string)-j+1)]

    # find unique substrings with number of repeats
    unique_substrings = {}
    for substring in substrings:
        if substring in unique_substrings:
            unique_substrings[substring] += 1
        else:
            unique_substrings[substring] = 1

    # only retain those with more than 3 repeats
    unique_substrings = {k: v for k, v in unique_substrings.items() if v > 3}

    # remove those with unbalanced brackets and not ending with )
    unique_substrings = {k: v for k, v in unique_substrings.items() if k.count('(') == k.count(')') and k.endswith(')')}

    # # ask for user input to remove unwanted substrings
    # index = input('Enter index of substring to remove (separated by comma): ')
    # index = [int(i) for i in index.split(',')]
    # for i in index:
    #     del unique_substrings[list(unique_substrings.keys())[i]]

    # # print the remaining substrings
    # print(unique_substrings)

    return unique_substrings


def find_pow(string, max_layers=1, min_repeats=2):
    '''
    Function to find all substrings of the form pow( ... ) in a string and return them as a dictionary with number of repeats.
    Grouped in terms of layers of pow( ... ) (i.e. pow( pow( ... ) ) is layer 2).
    max_layers is the maximum number of layers to consider.
    '''

    # find all substrings between pow( and ), making sure brackets are balanced
    substrings = []
    for i in range(len(string)):
        if string[i:i+4] == 'pow(':
            count = 1
            j = i+4
            while count > 0:
                if string[j] == '(':
                    count += 1
                elif string[j] == ')':
                    count -= 1
                j += 1
            substrings.append(string[i:j])

    # retain uniques and count number of repeats
    unique_substrings = {}
    for substring in substrings:
        if substring in unique_substrings:
            unique_substrings[substring] += 1
        else:
            unique_substrings[substring] = 1

    # retain those with at least min_repeats
    unique_substrings = {k: v for k, v in unique_substrings.items() if v >= min_repeats}
    
    # group by layer in a dictionary, where key is layer and value is list of substrings
    grouped = {}
    for substring in unique_substrings:
        # count number of layers in substring
        count = substring.count('pow(')
        if count in grouped:
            grouped[count].append(substring)
        else:
            grouped[count] = [substring]

    # only retain those with at most max_layers
    grouped = {k: v for k, v in grouped.items() if k <= max_layers}

    return grouped


def optimisation(string, max_layers=1, min_repeats=2):
    '''
    Function to replace c++ expressions with repeated pow() with a variable by using find_pow iteratively until no more repeated pow() are found.
    '''

    i = 1 # layer counter
    definition_string = '' # string to store variable definitions

    while True:
        # find all substrings of the form pow( ... ) in a string and group them by layer
        substrings = find_pow(string, max_layers=max_layers, min_repeats=min_repeats)

        print('Iteration {}:'.format(i))
        print(substrings)

        # check non-empty substrings
        if substrings == {}:
            break
        else:
            substrings = substrings[1]

        j = 1 # variable counter
        for substring in substrings:
            # replace substring with variable
            string = string.replace(substring, 'v{}{}'.format(i, j))

            # add variable definition to definition_string
            definition_string += 'long double v{}{} = {};\n'.format(i, j, substring)
            
            j += 1
        
        i += 1

    results = {'string': string, 'definition_string': definition_string, 'combined_string': definition_string + string}

    return results

TEST_STRING = r'''
    // gtt
	gtt = r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2));
	// grr
    grr = r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow(a52 + pow(r,2),-1)*pow((-2 + r)*r + pow(spin,2),-1);
	// gthth
    gthth = pow(r,2) + pow(spin,2)*pow(cos(th),2);
	// gpp
    gpp = pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2));
	// gtp
    gtp = -2*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2);
'''

# Test optimisation
print(find_repeats(TEST_STRING, min_length=4, max_length=30))
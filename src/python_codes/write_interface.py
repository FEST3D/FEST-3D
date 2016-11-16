import sys

def extract_method_name(codeline, method_type):
    '''Extract the method name from the given codeline.

    The method could be a subroutine or a function. The method name is present 
    immediately after the keyword 'subroutine' or 'function' and ends just 
    before the parenthesis containing the arguments begins.

    The input argument method_type should be either 'subroutine' or 'function'.
    '''
    if codeline.startswith(method_type):
        # Remove the method_type from the beginning
        # The + 1 accounts for the space immediately preceeding the name
        start = len(method_type) + 1
        # Find the parenthesis
        end = codeline.find('(')
        method_name = codeline[start:end]
    else:
        method_name = None
    return method_name

def extract_arg_types(code, line, arg_list, result_var):
    '''Extract the subroutine argument types.
    
    The lines following the method name and argument list will contain the 
    argument type declarations. This function will extract the argument types 
    and bundle them in a list. This list will contain lines that can be written 
    to the output file without modification.

    FORTRAN functions may also have a result variable. We need to get the type
    of this variable as well.

    This function returns a tuple containing the argument type list as its first
    item and the line number at which the last argument was found as its second 
    item.

    The function requires three inputs:
      code: The list of code lines
      line: The line number at which the subroutine name was found
      arg_list: The list of arguments to the subroutine
    '''
    arg_types = []
    if result_var:
        len_arg_list = len(arg_list) + 1
    else:
        len_arg_list = len(arg_list)
    # skip indicates whether the current line needs to be skipped or not
    skip = False
    # We are done when we obtain one type definition per argument in the list
    while len(arg_types) < len_arg_list:
        # Check if the current line is a variable declaration line
        # Each variable declaration line in fortran is of the form:
        #   type-specifier :: variable-list
        # First remove any trailing comments
        codeline = code[line].split('!')[0].strip()
        if codeline == 'interface':
            # The lines following this until the end interface statement is
            # reached need to be skipped.
            skip = True
        elif codeline == 'end interface':
            # Resume reading the lines (don't skip any more)
            skip = False
        if skip:
            line += 1
            continue
        if codeline.find('::') != -1:
            # This line is a variable declaration line
            # Separate the type-specifier from the variable list
            declaration_line = codeline.split('::')
            type_specifier = declaration_line[0].strip()
            var_list = declaration_line[1].strip()
            # Extract the variable names
            var_list = var_list.split(',')
            var_list = [var.strip() for var in var_list]
            # Check if the variable is a subroutine argument
            for var in var_list:
                if var in arg_list or var == result_var:
                    # This variable is an argument
                    # Save its type
                    type_info = type_specifier + ' :: ' + var
                    arg_types.append(type_info)
        line += 1
    return (arg_types, line)

def extract_method_args(code, line):
    '''Extract the method arguments.

    The arguments passed to the subroutine will be contained within parenthesis
    right after the subroutine name. These may also span several lines (lines in
    FORTRAN may be continued using the continuation character ('&')). 

    This function extracts all the arguments to a subroutine and bundles them in
    a list. The returned value is a tuple containing this list as its first item
    and the line number at which the last argument was found as its second item.

    The above is valid for functions as well.

    The function requires two inputs:
      code: The list of code lines
      line: The line number at which the subroutine name was found
    '''
    # Find the opening parenthesis
    arg_begin_index = code[line].find('(') + 1
    # Strip the subroutine name from the line
    arg_list = code[line][arg_begin_index:]
    # Split at the commas and strip spaces to get the arg names
    arg_list = arg_list.split(',')
    arg_list = [arg.strip() for arg in arg_list]
    # Check if the line is continued on the next line
    while arg_list[-1] == '&':
        # If it is, get the rest of the arguments from the next line
        line += 1
        temp_arg_list = code[line].split(',')
        temp_arg_list = [arg.strip() for arg in temp_arg_list]
        # Remove the '&' before extending the list
        arg_list.pop(-1)
        arg_list.extend(temp_arg_list)
    else:
        # If the current line does not continue to the next, there will be a
        # closing parenthesis at the end of the last argument. In the case of 
        # a function, there will be an additional result variable.
        # Remove these.
        end = arg_list[-1].find(')')
        # end now contains the index of the first closing parenthesis.
        # This is what we want
        arg_list[-1] = arg_list[-1][:end]
    return (arg_list, line)

def write_data(method_name, method_type, arg_list, result_var, arg_types):
    '''Write the interface file using the information obtained.'''
    # Create a buffer containing the output
    buf = method_type + ' ' + method_name
    buf += '('
    for arg in arg_list:
        buf += arg + ', '
    # Remove the final comma
    if buf[-2:] == ', ':
        buf = buf[:-2]
    buf += ')'
    if result_var and result_var != method_name:
        buf += ' result(' + result_var + ')'
    buf += '\n'
    for type_info in arg_types:
        buf += '    ' + type_info + '\n'
    buf += 'end ' + method_type + ' ' + method_name
    # Open the output file (to contain the interface)
    # Opening in write ('w') mode will truncate the file
    # This is fine since we will regenerate the file anyway
    output_filename = method_name + '.interface'
    output = open(output_filename, 'w')
    output.write(buf)
    output.close()
    print 'Interface written to: ' + output_filename

def extract_result_variable(codeline):
    '''Extract the function's result variable.

    Functions in FORTRAN can return one result. The variable for this is 
    declared using the result keyword. This function will extract the result
    variable name and return it.

    If the result variable is not given, the function name is used as the 
    result variable.
    '''
    start = codeline.find(' result(')
    if start == -1:
        print 'No result variable found.'
        return None
    # The start index of the variable is right after the ' result('.
    start += len(' result(')
    # The variable ends just before the closing parenthesis (after the 
    # corresponding opening parenthesis.
    end = codeline[start:].find(')') + start
    result_var = codeline[start:end]
    return result_var

def write_interface(codefile):
    '''Create a file with the FORTRAN interface for the subroutine / function.

    Updating a subroutine's arguments requires updating the interfaces in all 
    the files that call the subroutine. This can be avoided if these interfaces 
    are placed in a different file and included in the necessary FORTRAN code 
    files. This function can be used to create such interface files.

    The input codefile should comply with the usual styleguide; atleast the 
    following guidelines must be met:
      - The file should contain one or more subroutines.
      - Each subroutine should begin at the start of the line.
      - The arguments to the subroutine can be split across many lines. If this
        happens, the line continuation character (&) should be the last 
        character on any line that is to be continued.

    All the above is valid for functions as well.
    '''
    # Read the codefile
    try:
        f = open(codefile, 'r')
    except IOError:
        print 'File ' + codefile + ' not found.'
        return
    else:
        code = f.read()
        f.close()
    code = code.splitlines()
    # Traverse through the file looking for subroutines and functions
    line = 0
    while line < len(code):
        if code[line].startswith('subroutine'):
            method_type = 'subroutine'
        elif code[line].startswith('function'):
            method_type = 'function'
        else:
            line += 1
            continue
        # Extract the method name
        method_name = extract_method_name(code[line], method_type)
        print 'Found ' + method_type + ': ' + method_name
        # Get the method arguments
        (arg_list, line) = extract_method_args(code, line)
        # If the method is a function, get the result variable
        if method_type == 'function':
            result_var = extract_result_variable(code[line])
            if result_var is None:
                result_var = method_name
        else:
            result_var = None
        # Get the argument types
        (arg_types, line) = extract_arg_types(code, line, arg_list, result_var)
        # Write the interface to the output file
        write_data(method_name, method_type, arg_list, result_var, arg_types)
        # Proceed to the next line
        line += 1

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print 'Please specify atleast one file.'
        print 'Usage: python write_interface.py file-list'
    else:
        for arg in sys.argv[1:]:
            write_interface(arg)

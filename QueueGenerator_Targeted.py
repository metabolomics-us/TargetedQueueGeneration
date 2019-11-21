from pandas import DataFrame
from pandas import read_excel
from math import ceil
from math import floor
from random import shuffle
from os import path


def main():
    inputs = get_inputs()
    inj_info, max_sample, plate_num = make_samples(inputs)
    queue = assign_injection_names(sample_queue(max_sample, inputs, inj_info))
    queue_csv = DataFrame(data=queue, dtype=str)
    save_to_file(inputs, queue_csv)
    input("Press enter to exit program")


def assign_injection_names(queue: list):
    """
    After sample randomization and adding pools, generates dictionary as basis for dataframe
    :param queue: output of sample_queue()
    :return: dict consisting of lists for creating dataframe
    """
    i = 1
    final_queue = {'Sample Name': [],
                   'Well': [],
                   'Plate': [],
                   'Injection Name': [],
                   'InjWell': [],
                   'Non-random index': []}
    for sample in queue:
        sample_name = sample[0]
        final_queue['Sample Name'].append(sample_name)
        well = sample[1]
        final_queue['Well'].append(well)
        plate = sample[2]
        final_queue['Plate'].append(plate)
        final_queue['Injection Name'].append("{}-{}.{}-{}".format(i, well, plate, sample_name))
        final_queue['InjWell'].append(make_address(well, plate))
        final_queue['Non-random index'].append(sample[3])
        i += 1
    return final_queue


def make_address(well: str, plate: str):
    """
    Generate Vanquish Well ID from the plate ID and Well
    :param well: str, eg "C10"
    :param plate: str, If standard plate, will assign Blue plate, else use Green, Red, Yellow
    :return: str
    """
    plate_list = ['G', 'R', 'Y']
    try:
        plate = int(plate) % 3
        address = plate_list[plate] + ':' + well
        return address
    except ValueError:
        return 'B:' + well


def get_inputs():
    """ Returns dictionary with the following:
        List of samples
        Client name
        MX number
        Sample matrix
        Platform
        Random flag
        """
    inputs = dict()
    inputs['filepath'] = input("Please enter sample list filepath: ").strip().strip('"')
    df = read_excel(inputs['filepath'], header=None)
    inputs['client'] = df.at[2, 0][18:].split()[1]
    samples = list()
    inputs['samples'] = samples
    row_num, col_num = df.shape
    inputs['batch'] = get_batch_size()

    for ind in range(10, row_num):
        samples.append(df.at[ind, col_num-1].replace(" ", ""))
    inputs['minix'] = input("Please enter MX id: ").strip()
    inputs['platform'] = get_platform()
    inputs['matrix'] = df.at[8, 0].replace(
        "Specimen Type: (e.g.plasma, serum, stool…)",
        "").strip()
    inputs['random'] = get_bool('Do you want to randomize samples? (y/n): ')

    return inputs


def get_bool(prompt: str):
    """
    transforms a yes/no or y/n user input into a boolean
    :param prompt: text asking for
    :return:
    """
    while True:
        rand = input(prompt).strip().lower()
        if rand == 'y' or rand == 'yes':
            return True
        elif rand == 'n' or rand == 'no':
            return False
        else:
            print('Input invalid')


def get_platform():
    while True:
        try:
            print("Please indicate the study's platform:\n(1) Bile Acids\n(2) Steroids\n(3) Oxylipins")
            indicator = int(input("").strip()) - 1
            return ['BA', 'Ster', 'Oxy'][indicator]
        except (ValueError, IndexError):
            print('Invalid input. Please try again.')


def get_batch_size():
    """
    Fetches size of batches for use in generate_partition().
    Leaving blank returns 0
    returns int
    """
    while True:
        try:
            print('Please enter the number of samples between standards desired.')
            batch = input('(Leave blank for auto-calculation): ').strip()
            return int(batch)
        except ValueError:
            if batch == '':
                return 0
            else:
                print('Invalid response. Please enter an integer.')


def gen_wells():
    """
    Generates tuple with well names starting at A9 - H12
    :return: tuple of strings
    """
    letters = list('ABCDEFGH')
    wells = list()
    for row in letters:
        for col in range(1, 13):
            wells.append(row + str(col))
    return tuple(wells[8:])


def make_samples(inputs: dict):
    """
    :param inputs: dictionary output of get_inputs()
    :return:
        inj_info: list of tuples consisting of sample name, well, and plate they will end up in
        max_sample: int, length of samples
        plate_num: integer, # of plates needed for samples
    """
    samples = inputs['samples']
    well_list = gen_wells()
    max_sample = len(samples)
    inj_info = list()
    plate = list()
    plate_num = 1
    index = 0
    num = 1
    for sample in samples:
        plate.append(tuple([sample, well_list[index], plate_num, num]))
        index += 1
        num += 1
        if index >= 88:
            plate_num += 1
            index = 0
            if inputs['random']:
                shuffle(plate)
            inj_info.extend(plate)
            plate = list()
    if inputs['random']:
        shuffle(plate)
    inj_info.extend(plate)
    return inj_info, max_sample, plate_num


def find_cal_curve_number(sample_no: int):
    """
    :param sample_no: number of samples in the study
    :return: int, the number of cal_curves required from 0-8
    """
    x = sample_no + 25
    return ceil(x / 75)


def make_sample_partition(sample_no: int, batch_size: int):
    """
    :param sample_no: Number of samples in study
    :param batch_size: desired # of samples between cal curve instances
    :return: partition: tuple of int, distribution of the samples between cal curve instances
        If batch size was not specified, batch sizes will be roughly even (20-25 samples)
        If batch size was specified as n, the last batch size may from 1-n samples long
    """
    # Auto-calculation in case of no inputs
    if not batch_size:
        cal_curve_number = find_cal_curve_number(sample_no)
        batches = cal_curve_number * 3 - 1
        min_batch = floor(sample_no / batches)
        partition = [min_batch] * batches
        i = 0
        while i < sample_no % batches:
            partition[i] = min_batch + 1
            i += 1
        return tuple(partition)
    # Splitting into set batch size
    else:
        batches = floor(sample_no / batch_size)
        last_batch = sample_no % batch_size
        partition = [batch_size] * batches
        partition.append(last_batch)
        return tuple(partition)


def sample_queue(sample_no: int, inputs: dict, inj_info: list):
    """
    :param sample_no: int, number of client samples
    :param inputs: dict, user inputs, primarily needed for batch size and platform
    :param inj_info: list of sample tuples
    :return: list of tuples, represents the queue of injections including blanks, QCs, Stds, and samples
    """
    # Setting up parameters
    partition = make_sample_partition(sample_no, inputs['batch'])
    batch = 0
    batch_position = 0
    cal_curve_state = 2
    queue = list()
    platform = inputs['platform']
    cal_curve_no = 1
    pool_no = 1
    methanol_no = 1
    plate_num = 1
    # Begin with methanols, Std curves, and test injections
    queue.extend(methanols(3, methanol_no))
    methanol_no = 4
    queue.extend(cal_curve(0, 1, platform, methanol_no, pool_no))
    queue.extend(qc_blanks(plate_num))
    # Iterate variables
    methanol_no += 1
    pool_no += 1

    for sample in inj_info:
        # Check for adding in blanks and internal std checks at beginning of plate
        if type(sample[2]) is int and sample[2] != plate_num:
            plate_num += 1
            queue.extend(qc_blanks(plate_num))
        # Add sample
        queue.append(sample)
        # At the end of the batch, add cal curves & MeOH
        batch_position += 1
        if batch_position == partition[batch]:
            batch_position = 0
            batch += 1
            queue.extend(cal_curve(cal_curve_state, cal_curve_no, platform, methanol_no, pool_no))
            methanol_no += 1
            # Add on pool if necessary
            if cal_curve_state < 2:
                pool_no += 1
            # Iterating through cal curve
            if cal_curve_state == 3:
                cal_curve_state = 1
                cal_curve_no += 1
            else:
                cal_curve_state += 1
    # Final MeOH wash injections
    queue.extend(methanols(2, methanol_no))
    return queue


def qc_blanks(plate_num):
    blank_1 = ('Blank{}'.format(str(plate_num*2 - 1)), 'A1', str(plate_num), 'QC')
    blank_2 = ('Blank{}'.format(str(plate_num*2)), 'A2', str(plate_num), 'QC')
    ist1 = ('IST{}'.format(str(plate_num*2 - 1)), 'A3', str(plate_num), 'QC')
    ist2 = ('IST2'.format(str(plate_num*2)), 'A3', str(plate_num), 'QC')
    return [blank_1, blank_2, ist1, ist2]


def cal_curve(state: int, number: int, platform: str, methanol_no: int, pool_no: int):
    """
    Generates list of cal curve, qc, and  blank samples to append in sample_queue()
    :param state: int, represents the 0369, 036, 147, or 258 series of cal curve injs
    :param number: int, The number of cal curves already constructed, indicates row that it exists in
    :param platform: str, (BA Ster or Oxy)
    :param methanol_no: int, he number of methanols already run
    :param pool_no: int, pools already run
    :return: list of tuples
    """
    cal_list = list()
    std_plate_name = "Std" + str(ceil(number/8))

    row = list('ABCDEFG')[(number-1) % 8]
    if state < 2:
        wells = [0, 3, 6]
        if state == 0:
            wells.append(9)
    elif state == 2:
        wells = [1, 4, 7]
    else:
        wells = [2, 5, 8]

    for well in wells:
        inj_name = "{}Std{}".format(platform, str(well))
        well_name = row + str(well + 1)
        cal_list.append(tuple([inj_name, well_name, std_plate_name, 'QC']))
    cal_list.append(_methanol(methanol_no))
    if state < 2:
        cal_list.append(add_pool(pool_no))

    return cal_list


def _methanol(methanol_no: int):
    """
    Generates a single methanol injection tuple
    :param methanol_no: int, number of methanol injection
    :return: tuple, injection
    """
    meoh_wells = ['A11', 'A12', 'B11', 'B12', 'C11', 'C12', 'D11', 'D12', 'E11', 'E12', 'F11', 'F12', 'G11', 'G12']
    name = "MeOH-{}".format(str(methanol_no))
    well = meoh_wells[methanol_no % 12]
    plate = 'Std' + str(ceil(methanol_no/12))
    return tuple([name, well, plate, 'QC'])


def methanols(repl: int, methanol_no: int):
    """
    Generates repl number of methanols
    :param repl:
    :param methanol_no:
    :return: list of methanol tuples
    """
    methanols_list = list()
    for i in range(repl):
        methanols_list.append(_methanol(methanol_no + i))
    return methanols_list


def add_pool(pool_no: int):
    """
    Generates pool injection tuple
    :param pool_no: the number of pools injected
    :return: tuple
    """
    pool_wells = ['A5', 'A6', 'A7', 'A8']
    plasma = 1
    if plasma:
        name = 'Utak+IST{}'.format(pool_no)
    else:
        name = 'Pool{}'.format(pool_no)
    well = pool_wells[(pool_no-1) % 4]
    plate = ceil(pool_no / 4)
    return tuple([name, well, plate, 'QC'])


def save_to_file(inputs, queue_df):
    """
    Take the file path of the folder and save the new to_csv into the directory the original file is in.
    :return: None
    """
    save_file_name = "MX{}_{}_{}_Queue.csv".format(inputs['minix'], inputs['client'], inputs['platform'])
    directory = path.dirname(inputs['filepath'])
    queue_df.to_csv(path_or_buf=directory + "\\" + save_file_name)
    print('Queue saved to', directory + "\\" + save_file_name)


main()

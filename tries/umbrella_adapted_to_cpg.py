import numpy as np


def forward(observations):
    initial_probabilities = np.array([0.5, 0.5])
    forward_list = [initial_probabilities]
    rain_t_minus_1_given_observation_t_minus_1 = initial_probabilities

    # Looping through days
    index = 1
    for observation in enumerate(observations):

        # Calculating P(Xt | Et-1) by dotting P(Xt | Xt-1) with (Xt-1 | Et-1)
        rain_t_given_umbrella_minus_1 = np.dot(transition_probabilities, rain_t_minus_1_given_observation_t_minus_1)
        # print("Rain(X{}|e{}) = {}".format(index, (index - 1), rain_t_given_umbrella_minus_1))

        # If umbrella is observed, dot P(Et | Xt) with P(Xt | Et-1)
        if observation[1] == possible_observations['+']:
            rain_t_given_observation_t = np.dot(emission_probabilities_umbrella, rain_t_given_umbrella_minus_1)
        # If umbrella is not observed, dot P(~Et | Xt) with P(Xt | Et-1)
        else:
            rain_t_given_observation_t = np.dot(emission_probabilities_no_umbrella, rain_t_given_umbrella_minus_1)

        # Normalizing result, which is P(Xt | Et)
        rain_t_given_observation_t = rain_t_given_observation_t / rain_t_given_observation_t.sum()
        print("P(X{}|e1:{}) = {}".format(index, index, rain_t_given_observation_t))

        # Add to forward list for further processing in smoothing
        forward_list.append(rain_t_given_observation_t)

        # Set P(Xt-1 | Et-1) to P(Xt | Et). In other words, set current yesterday to today before next iteration
        rain_t_minus_1_given_observation_t_minus_1 = rain_t_given_observation_t

        index += 1

    return forward_list


def backward(observations, print_output):
    # Reversing observations
    observations = observations[::-1]

    # Initial value
    b_hat = np.array([1.0, 1.0])

    if print_output:
        print("Backwards {}:{} = {}".format(len(observations), len(observations), b_hat))

    # Set list of all backward calculations - P(Ek+1:t | Xk)
    backwards_list = [b_hat]

    # Looping through days
    for i, observation in enumerate(observations):

        # If umbrella is observed, dot P(Et | Xt) with P(Ek+1:t | Xk)
        if observation == possible_observations['+']:
            a = np.dot(emission_probabilities_umbrella, b_hat)
        # If umbrella is not observed, dot P(~Et | Xt) with P(Ek+1:t | Xk)
        else:
            a = np.dot(emission_probabilities_no_umbrella, b_hat)

        # Dotting a with P(Xt | Xt-1)
        b = np.dot(a, transition_probabilities)

        if print_output:
            print("Backwards {}:{} = {}".format(len(observations) - i - 1, len(observations), b))

        # Normalizing
        b = b / b.sum()

        # Adding to list for further processing in smoothing
        backwards_list.append(b)

        # Set last value to current value
        b_hat = b

    return backwards_list


def smoothing(forward_list, backward_list):
    # Reversing backwards list
    backward_list = backward_list[::-1]

    # Looping through backwards list
    for index, b in enumerate(backward_list):
        # Multiplying the element form the forward list with the element in the reversed backwards list
        rain_k_given_observation_e_1_t = np.multiply(forward_list[index], b)

        # Normalizing
        rain_k_given_observation_e_1_t = rain_k_given_observation_e_1_t / rain_k_given_observation_e_1_t.sum()

        if index > 0:
            print("P(X{}|e{}:{}) = {}".format(len(backward_list) - index, 1, len(forward_list) - 1,
                                              rain_k_given_observation_e_1_t))


possible_observations = {
    '+': True,
    '-': False
}

transition_probabilities = np.array([[0.9, 0.1],
                                     [0.05, 0.95]])

emission_probabilities_umbrella = np.array([[0.1477131, 0, 0, 0],
                                            [0, 0.3339165, 0, 0],
                                            [0, 0, 0.3646588, 0],
                                            [0, 0, 0, 0.1537116]])

emission_probabilities_no_umbrella = np.array([[0.2617500, 0, 0, 0],
                                               [0, 0.2470000, 0, 0],
                                               [0, 0, 0.2382500, 0],
                                               [0, 0, 0, 0.2530000]])

observations_1 = [
    "A", "G", "T", "G", "C", "C", "C", "A", "C", "C", "T", "C", "G", "A", "T", "G", "G", "G", "C", "C", "A", "A", "G", "A", "A", "T", "T", "G", "T", "G", "C", "G", "T", "C",
    "A", "G", "A", "A", "G", "C", "T", "T", "A", "G", "C", "A", "G", "C", "G", "C", "G", "G", "C", "T", "G", "T", "C", "G", "G", "C", "C", "A", "C", "C", "G", "A", "T", "G",
    "C", "A", "A", "A", "A", "T", "A", "C", "A", "T", "T", "A", "C", "A", "T", "C", "G", "C", "A", "A", "A", "T", "C", "A", "A", "G", "G", "C", "T", "C", "C", "C", "C", "A",
    "A", "T", "G", "G", "T", "T", "C", "G", "C", "T", "T", "C", "C", "C", "G", "G", "C", "T", "G", "C", "C", "C", "T", "C", "A", "G", "G", "T", "A", "C", "C", "G", "G", "C",
    "G", "A", "G", "C", "A", "A", "T", "G", "C", "G", "A", "G", "T", "A", "C", "T", "G", "T", "C", "G", "A", "A", "A", "C", "G", "C", "C", "A", "C", "A", "T", "T", "A", "A",
    "C", "T", "T", "T", "T", "A", "G", "C", "G", "G", "G", "G", "T", "A", "C", "C", "G", "G", "A", "T", "A", "G", "A", "A", "A", "A", "C", "A", "C", "C", "G", "C", "G", "T",
    "A", "T", "G", "C", "G", "G", "C", "C", "G", "G", "C", "T", "C", "G", "C", "A", "C", "G", "C", "A", "T", "G", "T", "A", "A", "A", "T", "A", "T", "A", "G", "T", "G", "A",
    "C", "C", "G", "C", "G", "G", "C", "C", "G", "G", "C", "A", "T", "A", "T", "C", "A", "G", "A", "C", "C", "C", "C", "A", "T", "G", "T", "A", "G", "A", "G", "G", "A", "G",
    "T", "A", "G", "T", "T", "T", "T", "T", "T", "A", "A", "G", "G", "G", "C", "G", "G", "G", "A", "T", "G", "C", "G", "G", "T", "C", "G", "G"
]


if __name__ == '__main__':

    print("Forward probabilities:")
    forward_1 = forward(observations_1)

    print("\nBackward probabilities")
    backward_1 = backward(observations_1, True)

    print("\nSmoothing probabilities")
    smoothing(forward_1, backward_1)

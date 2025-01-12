#plotting error
import matplotlib.pyplot as plt
import numpy as np
import pickle
import Magnus_Expansion as ME

def main():
    params = ['A','B','C']
    for param in params:
        with open("Data/{}_dynamics_closed_sixth_order".format(param) ,"rb") as file:
            sixth_order_ME_data = pickle.load(file)
        with open("Data/{}_dynamics_closed_fourth_order".format(param),"rb") as file:
            fourth_order_ME_data = pickle.load(file)
        with open("Data/{}_dynamics_closed_RK4".format(param),"rb") as file:
            fourth_order_RK_data = pickle.load(file)
        with open("Data/{}_dynamics_closed_RK5".format(param),"rb") as file:
            fifth_order_RK_data = pickle.load(file)
        print(type(fourth_order_ME_data[1][0][0]))
        print(fourth_order_ME_data[1][0][0].shape)
        print("ME error is",ME.error_matrix(sixth_order_ME_data[1][0][0][-1],fourth_order_ME_data[1][0][0][-1]))
        fifth_order_RK_data_np = convert_lists_to_array(fifth_order_RK_data[1])
        #print(fifth_order_RK_data_np[0])
        print("RK4 error is", ME.error_matrix(sixth_order_ME_data[1][0][0][-1],fourth_order_RK_data[1][-1]))
        print("RK5 error is", ME.error_matrix(sixth_order_ME_data[1][0][0][-1],fifth_order_RK_data_np[-1]))
        print("RK4 trace is", get_traces(fourth_order_RK_data[1][-1]))
        print("RK trace is ",get_traces(fifth_order_RK_data_np[-1]))
        print("ME_4th trace is ", get_traces(fourth_order_ME_data[1][0][0][-1]))
        print("ME_6th trace is ", get_traces(sixth_order_ME_data[1][0][0][-1]))
        
    for param in params:
        with open("Data/{}_dynamics_sixth_order".format(param) ,"rb") as file:
            sixth_order_ME_data = pickle.load(file)
        with open("Data/{}_dynamics_fourth_order".format(param),"rb") as file:
            fourth_order_ME_data = pickle.load(file)
        with open("Data/{}_dynamics_RK4".format(param),"rb") as file:
            fourth_order_RK_data = pickle.load(file)
        with open("Data/{}_dynamics_RK5".format(param),"rb") as file:
            fifth_order_RK_data = pickle.load(file)
        print(type(fourth_order_ME_data[1][0][0]))
        print(fourth_order_ME_data[1][0][0].shape)
        print(ME.error_matrix(fourth_order_ME_data[1][0][0][-1],sixth_order_ME_data[1][0][0][-1]))
        fifth_order_RK_data_np = convert_lists_to_array(fifth_order_RK_data[1])
        #print(fifth_order_RK_data_np[0])
        print(ME.error_matrix(fifth_order_RK_data_np[-1],sixth_order_ME_data[1][0][0][-1]))
        print("RK trace is ",get_traces(fifth_order_RK_data_np[-1]))
        print("ME_4th trace is ", get_traces(fourth_order_ME_data[1][0][0][-1]))
        print("ME_6th trace is ", get_traces(sixth_order_ME_data[1][0][0][-1]))
        
def get_traces(unitary_operator):
    #finds the trace for 3 initial conditions
    traces = np.zeros(3,dtype=np.complex128)
    for i in range(3):
        traces[i] = np.trace(unitary_operator[:,i*4].reshape(3,3)) 
    return traces
    
def convert_lists_to_array(lists_data):
    #takes data stored in lists of lists and converts them to numpy array
    #used to deal with data produced by scipy.RK45
    array_data = np.array(lists_data)
    print(array_data.shape)
    return array_data.T.reshape(8193,9,9)
if __name__ == "__main__":
    main()
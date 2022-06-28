# -*- coding: utf-8 -*-
"""
Created on Wed May  4 11:17:19 2022

@author: tbart
"""

# Import the .NET class library
import clr

# Import python sys module
import sys

# Import os module
import os

# Import System.IO for saving and opening files
from System.IO import *

# Import C compatible List and String
from System import String
from System.Collections.Generic import *
from System.Threading import AutoResetEvent



# Add needed dll references
sys.path.append(os.environ['LIGHTFIELD_ROOT'])
sys.path.append(os.environ['LIGHTFIELD_ROOT']+"\\AddInViews")
clr.AddReference('PrincetonInstruments.LightFieldViewV5')
clr.AddReference('PrincetonInstruments.LightField.AutomationV5')
clr.AddReference('PrincetonInstruments.LightFieldAddInSupportServices')

# PI imports
from PrincetonInstruments.LightField.Automation import Automation
from PrincetonInstruments.LightField.AddIns import CameraSettings
from PrincetonInstruments.LightField.AddIns import DeviceType
from PrincetonInstruments.LightField.AddIns import ExperimentSettings

def convert_buffer(net_array, image_format):
    src_hndl = GCHandle.Alloc(net_array, GCHandleType.Pinned)
    try:
        src_ptr = src_hndl.AddrOfPinnedObject().ToInt64()

        # Possible data types returned from acquisition
        if (image_format==ImageDataFormat.MonochromeUnsigned16):
            buf_type = ctypes.c_ushort*len(net_array)
        elif (image_format==ImageDataFormat.MonochromeUnsigned32):
            buf_type = ctypes.c_uint*len(net_array)
        elif (image_format==ImageDataFormat.MonochromeFloating32):
            buf_type = ctypes.c_float*len(net_array)
                    
        cbuf = buf_type.from_address(src_ptr)
        resultArray = np.frombuffer(cbuf, dtype=cbuf._type_)

    # Free the handle 
    finally:        
        if src_hndl.IsAllocated: src_hndl.Free()        
    # Make a copy of the buffer
    return np.copy(resultArray)

def ManipulateImageData(dat, buff):
    im = convert_buffer(dat,buff.Format)
    im = np.reshape(im,(height,width))
    return im
    
def AcquireMoveAndLock():
    experiment.ImageDataSetReceived += experimentDataReady
    print("Acquiring...")  
    experiment.Acquire()            
    acquireCompleted.WaitOne()

def experimentDataReady(sender, event_args):
    print("image frames acquired")
    if (experiment.IsReadyToRun):
        global i
        frames = event_args.ImageDataSet.Frames
        i+=frames   #in case an event returns multiple frames
        buffer = event_args.ImageDataSet.GetFrame(0, frames-1) #1st ROI and most recent frame in the event
        image_data = buffer.GetData()
        data_array = ManipulateImageData(image_data,buffer)
        print("cool, did some stuff!")
        event_args.ImageDataSet.Dispose()

def set_value(setting, value):    
    # Check for existence before setting
    # gain, adc rate, or adc quality
    if experiment.Exists(setting):
        experiment.SetValue(setting, value)
        
def experiment_completed(sender, event_args):
    global image_data, data_array
    print("...Acquisition Complete! I think")   
    acquireCompleted.Set()
    experiment.ImageDataSetReceived -= experimentDataReady
    

def device_found():
    # Find connected device
    for device in experiment.ExperimentDevices:
        if (device.Type == DeviceType.Camera):
            return True
     
    # If connected device is not a camera inform the user
    print("Camera not found. Please add a camera and try again.")
    return False  


        
# Create the LightField Application (true for visible)
# The 2nd parameter forces LF to load with no experiment 
auto = Automation(True, List[String]())
experiment = auto.LightFieldApplication.Experiment
acquireCompleted = AutoResetEvent(False)


# Get experiment object
experiment.Load("pydemo")
experiment.ExperimentCompleted += experiment_completed
#experiment.SetValue(ExperimentSettings.AcquisitionFramesToStore, 11)



if (device_found()==True): 
    print("device found")
    #Set exposure time
    set_value(CameraSettings.ShutterTimingExposureTime, 20.0)

    AcquireMoveAndLock();
else: 
    print("nothing will happen")
    




""" 
This module implements a couple of custom exceptions. Their names are a bit more descriptive than
default exceptions.
"""
class MalformedDataError(Exception):
    "Raise this exception for data that are malformed"

class RequestError(Exception):
    "Raise this exception when an API call returns a non-200 code"

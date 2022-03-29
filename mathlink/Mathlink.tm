:Begin:
:Function:       set_alpha_precision
:Pattern:        SetAlphaPrecision[precision_Real]
:Arguments:      { precision }
:ArgumentTypes:  { Real }
:ReturnType:     Manual
:End:

:Begin:
:Function:       fit_model
:Pattern:        PowerlawFitModel[data_List, distributionType_String]
:Arguments:      { data, distributionType }
:ArgumentTypes:  { IntegerList, String }
:ReturnType:     Manual
:End:

:Begin:
:Function:       fit_model
:Pattern:        PowerlawFitModel[data_List, xParameter_Integer, distributionType_String]
:Arguments:      { data, xParameter, distributionType }
:ArgumentTypes:  { IntegerList, Integer, String }
:ReturnType:     Manual
:End:

:Begin:
:Function:       calculate_gof
:Pattern:        PowerlawGoodnessOfFit[data_List, replicas_Integer]
:Arguments:      { data, replicas }
:ArgumentTypes:  { IntegerList, Integer }
:ReturnType:     Manual
:End:

:Begin:
:Function:       calculate_gof
:Pattern:        PowerlawGoodnessOfFit[data_List, replicas_Integer, xParameter_Integer, distributionType_String, syntheticGeneratorMode_String]
:Arguments:      { data, replicas, xParameter, distributionType, syntheticGeneratorMode }
:ArgumentTypes:  { IntegerList, Integer, Integer, String, String }
:ReturnType:     Manual
:End:
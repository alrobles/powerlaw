:Begin:
:Function:       set_alpha_precision
:Pattern:        SetAlphaPrecision[precision_Real]
:Arguments:      { precision }
:ArgumentTypes:  { Real }
:ReturnType:     Manual
:End:

:Begin:
:Function:       set_test_statistic
:Pattern:        SetTestStatistic[testStatistic_String]
:Arguments:      { testStatistic }
:ArgumentTypes:  { String }
:ReturnType:     Manual
:End:

:Begin:
:Function:       fit_model
:Pattern:        PowerlawFitModel[data_List]
:Arguments:      { data }
:ArgumentTypes:  { IntegerList }
:ReturnType:     Manual
:End:

:Begin:
:Function:       fit_model
:Pattern:        PowerlawFitModel[data_List, xMin_Integer]
:Arguments:      { data, xMin }
:ArgumentTypes:  { IntegerList, Integer }
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
:Pattern:        PowerlawGoodnessOfFit[data_List, xMin_Integer, replicas_Integer]
:Arguments:      { data, xMin, replicas }
:ArgumentTypes:  { IntegerList, Integer, Integer }
:ReturnType:     Manual
:End:
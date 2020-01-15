# Post Processing of Fluent .dat Files with CPAD_oD function

Just change name ``mypath`` in file ``write_eval_scheme.py`` to current Fluent directory.
Then run the file with redirection of output to journal file.

```
python write_eval_scheme.py > cpad_pp.jou
```

Run the generated Fluent journal file, with the compiled cpad_udf_library, inside ANSYS Fluent.
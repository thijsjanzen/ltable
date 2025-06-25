# ltable

Ltables describe the procession of new species arising through speciation,
where each species is tracked in the table on a separate row, indicating the
time of birth, parent ID, self ID and finally the time of death. If the species
did not die (go extinct), the time of death is -1.

Ltables can be converted to phy objects using either `[DDD::L2phylo]` or 
`[treestats::l_to_phylo]`. However, by this conversion, some information is lost,
primarily the information on parentage. Instead, this R package allows for the 
direct visualisation of the Ltable, leaving information of parentage available.

## Example
```
res <- ltable::sim_bd(parameters = c(0.5, 0),
                      crown_age = 5,
                      seed = 4)
ltable::plot_ltable(res$ltable)
```
Yields:

<img src="https://github.com/user-attachments/assets/b7301104-2c03-4235-b0aa-201cb55c9bf2" alt="drawing" width="400"/>

Which represents the same phylogeny as:

<img src="https://github.com/user-attachments/assets/78752c8e-6346-432f-85f0-f5c25163ced2" width = "400"/>


# Writing Guidelines and Revision Plan

## Purpose
This document serves as a reference for polishing the methods section of the thesis. The introduction has already provided high-level explanations of time-varying DCM approaches and their applications. The methods section focuses on technical detail and implementation specifics.

## Style Guidelines

### Language and Spelling
- **British spelling throughout**: Use "modelling" (not "modeling"), "optimised" (not "optimized"), "parameterised" (not "parameterized"), "favour" (not "favor"), etc.
- **No em dashes**: Replace with semicolons, full stops, or restructure sentences for clarity
- **Minimal hedging**: Prefer direct technical statements over qualifiers like "arguably", "perhaps", "it could be said that"

### Tone and Precision
- **Technical, not verbose**: Write concisely with precise terminology
- **Direct statements**: Avoid unnecessary introductory phrases
- **Mathematical rigour**: Ensure notation is consistent and clearly defined
- **Biological grounding**: Technical claims about neural mechanisms should be specific and accurate

### Structural Principles
- **Paragraph length**: Break up paragraphs exceeding 150 words; each should have a clear focus
- **Subsection hierarchy**: Use clear heading structures to guide the reader
- **Transitions**: Ensure logical flow between sections without excessive signposting
- **Equations**: Introduce mathematical notation before use; explain variables clearly

## Section-by-Section Plan

1. Time varying DCM
1.1. Intro technical
    - fast and slow decomposition
1.2. Basis DCM
    - technical primer
1.3. Sequential DCM
What does seqDCM introduce?
- ability to model timevarying Parameter Trajectories of all neural parameters
    - biological interpretation of slow-varying parameters
- continous data modelling and interevent transitions accounoting for carryover effects and non-linear hystersis effects
- can model both exogenous and endogenous transitions
- irregular transitions
- complexity


2. Sequential DCM implementation
    2.1 Trajectory Function 
        2.1.1 mathematical implementation
        2.1.2 log-space parameterisation
        2.1.3 covaraince
    2.2 Inversion Scheme
        2.2.1 generic DCM extension
        2.2.2 incorpotin trajecotry function into inputs for inversion
        2.2.3 input functions scp
        2.2.4 parellisation
    2.3  Scaling comparisson between basisDCM and sequentialDCM




## Revision Process

For each section, the revision will:
1. **Identify specific issues**: Line-by-line review noting problems
2. **Propose improvements**: Specific edits with rationale
3. **Verify consistency**: Check cross-references, notation, terminology
4. **Confirm British spelling**: Systematic check
5. **Remove em dashes**: Replace with appropriate punctuation

## Success Criteria

A successfully revised section will:
- Use consistent British spelling
- Contain no em dashes
- Have focused paragraphs (typically 50-100 words)
- Present technical content precisely without excessive verbosity
- Flow logically from one idea to the next
- Define all mathematical notation clearly
- Support biological claims with appropriate specificity

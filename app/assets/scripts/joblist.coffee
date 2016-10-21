sortObjectsArray = (objectsArray, sortKey) ->
# Quick Sort:
  retVal = undefined
  if 1 < objectsArray.length
    pivotIndex = Math.floor((objectsArray.length - 1) / 2)
    # middle index
    pivotItem = objectsArray[pivotIndex]
    # value in the middle index
    less = []
    more = []
    objectsArray.splice pivotIndex, 1
    # remove the item in the pivot position
    objectsArray.forEach (value, index, array) ->
      if value[sortKey] <= pivotItem[sortKey] then less.push(value) else more.push(value)
      return
    retVal = sortObjectsArray(less, sortKey).concat([ pivotItem ], sortObjectsArray(more, sortKey))
  else
    retVal = objectsArray
  retVal


@a = ['0', 'p', 'q', 'r', 'e', 'd','i']
@jobs = {}


class Job
  constructor: (params) ->
    @mainID =  params.mainID
    @toolname = params.toolname
    @job_id = m.prop params.job_id
    @state = m.prop params.state


jobs.JobList = Array

jobs.vm = do ->
  vm = {}

  vm.init = ->

    m.request({url: "/api/jobs", method: "GET"}).then (jobs) ->
      vm.list = jobs.map (job) -> new Job(job)


    vm.stateList = {"0": "Partially Prepared", "p": "Prepared", "q": "Queued", "r": "Running", "e": "Error", "d": "Done", "i": "Submitted"}

    # Remove on Job with a certain mainID from the JObList
    vm.delete = (mainID, fromServer=false) ->

      oldLen = vm.list.length
      vm.list = vm.list.filter (job) -> job.mainID != mainID
      if vm.list.length < oldLen
        sendMessage("type": (if fromServer then "DeleteJob" else "ClearJob") ,"mainID":mainID)


    # Update a Job Object
    vm.update = (receivedJob) ->
      updatedJob = new Job(receivedJob)
      i = 0
      while i < vm.list.length
        job = vm.list[i]
        if job.mainID == updatedJob.mainID
          vm.list[i] = updatedJob
          return
        i++
      vm.list.push(updatedJob)


    # Sort by toolname // TODO the sorting algorithm works (tested it with the same JSON data in a plain html file), the sorting here does not stabilize, though
    vm.alphaSort = () ->
      vm.list = sortObjectsArray(vm.list, 'toolname')

    # Sort by job id
    vm.numericSort = () ->
      vm.list = sortObjectsArray(vm.list, 'job_id')


    vm.getJobState = (receivedJob) ->

      jsonString = JSON.stringify(receivedJob)
      if jsonString.indexOf('\"state\":3') > -1
        return 'running'
      else if jsonString.indexOf('\"state\":4') > -1
        return 'error'
      else if jsonString.indexOf('\"state\":5') > -1
        return 'done'
      else
        return 'other'

    vm.getLastJob = () ->
      return vm.list[vm.list.length-1]

  vm
#the controller defines what part of the model is relevant for the current page
#in our case, there's only one view-model that handles everything

jobs.controller = ->
  jobs.vm.init()

  select: (all) ->
    $('input:checkbox.sidebarCheck').each ->
      $(this).prop 'checked', if all then "checked" else ""

  delete: (fromServer) ->

    message = if fromServer then "Do you really want to delete the selected jobs permanently?"
    else "Do you really want to clear the selected jobs from the joblist?"
    ids = ($("input:checkbox.sidebarCheck:checked").map () -> $(this).val())
    if ids.length > 0 and confirm(message)
      for id in ids
        jobs.vm.delete(id, fromServer)

  toolsort: -> jobs.vm.alphaSort()

  idsort: -> jobs.vm.numericSort()



#here's the view
jobs.view = (ctrl) -> [

  m "div", {class: "button job-handle"}, [

    m "div", {class: "delete", onclick: ctrl.delete.bind(ctrl, true)}, "Delete"
    m "div", {onclick: ctrl.delete.bind(ctrl, false)} ,"Clear"
  ]

  m "div", {class: "button job-button"}, [

    m "div", [
      m "span",  {onclick: ctrl.select.bind(ctrl, true)}, "All"
      m "span", "/"
      m "span",  {onclick: ctrl.select.bind(ctrl, false)}, "None"
    ]
    m "div", {class: "idsort", onclick: ctrl.idsort}, "ID"
    m "div", {class: "toolsort", onclick: ctrl.toolsort}, "Tool"
  ]

  jobs.vm.list.map (job) ->
    m "div", {class: "job #{a[job.state()]}"},  [

      m "div", {class: "checkbox"}, [
            m "input", {type: "checkbox", class: 'sidebarCheck', name: job.mainID, value: job.mainID, id: job.mainID}
            m "label", {for: job.mainID}
        ]
      m "div", {class: "jobid"},  m 'a[href="/#/jobs/' + job.mainID + '"]', job.job_id()

      m "div", {class: "toolname"}, job.toolname.substr(0,4)
    ]
]


###
  <div class="button job-handle" style="display: flex;">
            <div style=" border-right: 1px solid lightgray;" onclick="deleteIDs();" data-tooltip title="Delete selected jobs from database">Delete</div>
            <div onclick="clearIDs();" data-tooltip title="Clear selected jobs from joblist">Clear</div>

        </div>
###



